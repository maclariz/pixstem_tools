import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt 
from matplotlib.ticker import MultipleLocator
from alive_progress import alive_bar

class stem:
  def discfloat(ci, cj, rmin, rmax, segments):
    # This is a helper function for polar transformation which generates a meshgrid
    # of float values corresponding to r and phi values back in a cartesian i, j frame.  
    # These can (after rounding to integers) then be used to lookup all the right pixels 
    # in a diffraction pattern to allow quick transformation from cartesian to polar
    # representations of the data using fancy lookup (i.e. passing lists of array indices
    # to slice an array).
    # I have chosen to use i and j to refer to vertical and horizontal axes in data
    # (rather than x and y, because some packages use x for horizontal,
    # and some use x for axis 0, which is vertical in python arrays)
    
    coords = np.zeros(shape=(2,rmax-rmin,segments))
    phi = np.arange(0, 2*np.pi, 2*np.pi/segments)
    r = np.arange(rmin,rmax)
    r_phi_mesh = np.meshgrid(r,phi)
    i = -r_phi_mesh[0]*np.sin(r_phi_mesh[1]) + ci
    j = r_phi_mesh[0]*np.cos(r_phi_mesh[1]) + cj
    return np.array([i, j])
    
    # returns array of dimensions:
      # 0 is the switch of i or j
      # 1 is the azimuth dimension
      # 2 is the radius / two theta dimension

  def polarttransform(DP, ci, cj, rmin, rmax, segments, simple=True):
    # Performs a polar transform of one single diffraction pattern
    disc = discfloat(ci, cj, rmin, rmax, segments) # get basic disc of all transform positions
    # Note: if segments are too few, then each segment will cover too many pixels in 
    # original dataset, simple mapping will be a very poor representation and even
    # simple == False mapping won't really represent all intensity in that patch of the original
    # You are better off using an appropriate number of segments to approximately match the 
    # 2*pi*r at the largest radius of interest in your analysis to get a good sampling of the 
    # original data in your transform
    
    if simple==True:
        # simply works things out with nearest position to the r, theta positions mapped
        # back into Cartesian
        pos = np.round(disc,0).astype('int16') # round the disc array to nearest integer
        pos2 = pos.reshape((2,pos.shape[1]*pos.shape[2])) # turn to a 1D list
        pt = DP[pos2[0],pos2[1]].reshape(pos.shape[1],pos.shape[2]).T
        # calculate PT by using the pos2 array to slice the original array

    if simple==False:
        # works out the weighted average of the four pixels surrounding the float 
        # r, theta positions mapped back into Cartesian
        shape = disc[0].shape[0]*disc[0].shape[1]
        disc0 = disc[0].reshape(shape) # turn into linear array of i positions
        disc1 = disc[1].reshape(shape) # turn into linear array of j positions
        ui = np.floor(disc0).astype('int16') # find upper i pixel array
        li = np.ceil(disc0).astype('int16') # find lower i pixel array
        li = np.where(li==ui,li+1,li) # deals with the case of an exact hit on an i position
        lj = np.floor(disc1).astype('int16') # find left j pixel array
        rj = np.ceil(disc1).astype('int16') # find right j pixel array
        rj = np.where(rj==lj, lj+1, rj) # deals with the case of an exact hit on a j position
        wul = (1-(disc0-ui))*(1-(disc1-lj)) # weighting parameter upper left
        wur = (1-(disc0-ui))*(1-(rj-disc1)) # weighting parameter upper right
        wll = (1-(li -disc0))*(1-(disc1-lj)) # weighting parameter lower left
        wlr = (1-(li -disc0))*(1-(rj-disc1)) # weighting parameter lower right
        pt = (
            DP[ui,lj]*wul +
            DP[ui,rj]*wur +
            DP[li,lj]*wll +
            DP[li,rj]*wlr
        ).reshape(disc[0].shape[0],disc[0].shape[1]).T
        
    # Now weight result by pixel area in transform image
    radweight = np.arange(rmin, rmax)*2*np.pi/segments
    azi = np.ones(shape=(segments))
    rweighting = np.meshgrid(azi,radweight)[1]
    PTDP = pt*rweighting
    
    return PTDP

  def PT4Dinone(dataset, ci, cj, rmin, rmax, segments, simple=True):
    # This function runs the polar transform over an entire 4DSTEM dataset
    # dataset is a 4DSTEM dataset as a numpy array
    # dimensions 0 and 1 are the vertical and horizontal dimensions of the image
    # dimensions 2 and 3 are the vertical and horizontal dimensions of the diffraction patterns
    # ci and cj are the pattern centres, either as floats or as arrays of floats to match dataset
    # rmin and rmax are the minimum and maximum radii for the output transform
    # you will save time and memory if you do not transform everything but focus on the radii of most interest
    # segments is the number of segments in azimuthal angle to split into (e.g. 360 gives 1 degree segments)
    # simple == True just calculates mapping one pixel in dataset to one in output
    # simple == False maps a weighted average of four pixels in dataset to one in output
    
    Ri, Rj = dataset.shape[0], dataset.shape[1]
    with alive_bar(Ri*Rj, force_tty=True) as bar:
      PT4D = np.zeros(shape=(Ri,Rj,rmax-rmin,segments))

      # version of calculation for a single value for pattern centre (probably okay if there is not much 
      # movement of pattern centre in dataset).
      
      if isinstance(ci, int):
            for i in range(Ri):
                for j in range(Rj):
                    PT4D[i,j,:,:] = polarttransform(
                        dataset[i,j,:,:], 
                        ci, cj, 
                        rmin, rmax, 
                        segments, 
                        simple=simple
                    )
                    bar()
            return PT4D
    
        # version of calculation for an array of pattern centres (more robust of there are some descan issues)
            
        elif isinstance(ci, np.ndarray):
            if ci.shape[0]==Ri and cj.shape[1]==Rj:
                for i in range(Ri):
                    for j in range(Rj):
                        PT4D[i,j,:,:] = polarttransform(
                            dataset[i,j,:,:], 
                            ci[i,j], cj[i,j], 
                            rmin, rmax, 
                            segments, 
                            simple=simple
                        )
                        bar()
                return PT4D
            else:
                print('The array size for the pattern centres does not match the dataset')
                return

  def plotpolar(polar, rmin, rmax, lines, title):
    # a little function for just plotting polar transformed datasets as a sanity check
    # no return, just an inline plot with appropriate angle labels
    
    # Parameters:
    # polar: a single polar transformed diffraction pattern
    # rmin: minimum radius of the transform in pixels
    # rmax: maximum radius of the transform in pixels
    # lines: a set of five line positions to delineate the Laue zone, between 2 and 3 is used
    # in this notebook for the fitted area
    fig, ax = plt.subplots(figsize=(12,6))
    ax.imshow(
        polar, 
        vmin=np.percentile(polar,5), 
        vmax=np.percentile(polar,98),
        cmap='turbo',
        extent = [0, 360, rmax, rmin],
        aspect=3
    ) #plotting transform 
                                                                  #with some background removed
    ax.set_xlabel(r'$\phi,$ deg')
    ax.set_ylabel('radius, pixels')
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.xaxis.set_major_locator(MultipleLocator(60))

    ax.hlines(lines+rmin,0,359,color=['m','k','w','w','m']) #adding lines to mark out position of HOLZ ring
    ax.set_title(title)

  def fun_cos_sq(phi,A2,phi2,A1,phi1,B):
  #intensity function for fitting of polar transformed data
    return (A2*np.cos(np.radians(phi-phi2))**2+
            A1*np.cos(np.radians(phi-phi1))+
            B)

def fitIntensitypattern(
    func, 
    pattern, 
    p0=[2000,90,2000,0,1000], 
    bounds=([0, 0, 0, -180, 0], [np.inf, 180, np.inf, 180, np.inf])
):
    ''' 
    fit HOLZ ring intensity to a sinusoidal function in a single pattern
    
    Parameters
    ----------
    
    func : sinusoidal function to fit data to
    
    pattern: 1D numpy array of intensity values as a function of i, j and azimuthal angle
    
    
    Returns
    -------
    
    fitParams : 3D numpy array of fit parameters (i_max x j_max x 5 array)
    
    fitCov: 4D numpy array with covariance matrix of fitted parameters (i_max x j_max x 5 x 5 array)

    Note
    ----

    p0 and bounds values are set to appropriate value for some experimental datasets of mine, but won't
    be suitable for everything.
    A2, A1 and B can range from 0 to infinity but a suitable starting value is best found by fitting one 
    pattern and checking these parameters before fitting a dataset.
    phi2 ranges in a 180 degree range (being a 2-fold function).  It's up to you what is the most sensible range.  
    If the value is significantly away from 0, then 0-180 is probably sensible.  But if around 0, then -90 to 90
    is probably better.
    phi1 ranges in a 360 degree range, and it's up to you to work out if 0 - 360 or -180 - 180 (or some other range)
    is more sensible.
    '''
    xDim, yDim, zDim = data.shape
    fitParams = np.zeros((xDim,yDim,5))
    fitCov = np.zeros((xDim,yDim,5,5))

    with alive_bar(yDim*xDim, force_tty=True) as bar:
        for i in range(yDim):
            for j in range(xDim):
                pop, pcov = op.curve_fit(
                    fun_cos_sq, 
                    np.arange(zDim), 
                    data[j,i], 
                    p0=p0, 
                    bounds=bounds
                )
                fitParams[j,i] = pop
                fitCov[j,i] = pcov
                bar()
    return fitParams, fitCov
  
def fitIntensitydataset(
    func, 
    data, 
    p0=[2000,90,2000,0,1000], 
    bounds=([0, 0, 0, -180, 0], [np.inf, 180, np.inf, 180, np.inf])
):
    ''' 
    fit HOLZ ring intensity to a sinusoidal function over entire 3D dataset
    
    Parameters
    ----------
    
    func : sinusoidal function to fit data to
    
    data: 3D numpy array of intensity values as a function of i, j and azimuthal angle
    
    
    Returns
    -------
    
    fitParams : 3D numpy array of fit parameters (i_max x j_max x 5 array)
    
    fitCov: 4D numpy array with covariance matrix of fitted parameters (i_max x j_max x 5 x 5 array)

    Note
    ----

    p0 and bounds values are set to appropriate value for some experimental datasets of mine, but won't
    be suitable for everything.
    A2, A1 and B can range from 0 to infinity but a suitable starting value is best found by fitting one 
    pattern and checking these parameters before fitting a dataset.
    phi2 ranges in a 180 degree range (being a 2-fold function).  It's up to you what is the most sensible range.  
    If the value is significantly away from 0, then 0-180 is probably sensible.  But if around 0, then -90 to 90
    is probably better.
    phi1 ranges in a 360 degree range, and it's up to you to work out if 0 - 360 or -180 - 180 (or some other range)
    is more sensible.
    '''
    xDim, yDim, zDim = data.shape
    fitParams = np.zeros((xDim,yDim,5))
    fitCov = np.zeros((xDim,yDim,5,5))

    with alive_bar(yDim*xDim, force_tty=True) as bar:
        for i in range(yDim):
            for j in range(xDim):
                pop, pcov = op.curve_fit(
                    fun_cos_sq, 
                    np.arange(zDim), 
                    data[j,i], 
                    p0=p0, 
                    bounds=bounds
                )
                fitParams[j,i] = pop
                fitCov[j,i] = pcov
                bar()
    return fitParams, fitCov
