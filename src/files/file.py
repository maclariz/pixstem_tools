import os, fnmatch
# gives you the list of files fitting some spec in your working directory
class files: 
	def file_list(path, searchterm = False):
  
	# path is a normal path, ending with a /
	# searchterm is False for no specific term, or a text string you want to find 
	# in the filenames to be listed, which will need to include wildcards
		files = []
    		paths = os.listdir(path)
    		if searchterm == False:
        		files = paths
        		return sorted(files)
    		elif isinstance(searchterm, str):
        		for file in paths:
            		if fnmatch.fnmatch(file, searchterm): 
                	#taking only files in dataset that match the search term
                	files += [file]
        		return sorted(files)
    		else:
        		print('Error: use a text string as search term!')
        		return
