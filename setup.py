from setuptools import setup, find_packages
from distutils.util import convert_path

with open("README.md", "r") as f:
    long_description = f.read()

version_ns = {}
vpath = convert_path("version.py")
with open(vpath) as version_file:
    exec(version_file.read(), version_ns)

setup(
    name="pixstem_tools",
    version=version_ns["__version__"],
    packages=find_packages(),
    description="An open source python package for 4DSTEM data manipulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/maclariz/pixstem_tools/",
    author="Ian MacLaren",
    author_email="ian.maclaren@glasgow.ac.uk",
    license="GNU GPLv3",
    python_requires=">=3.60",
    install_requires=[
        "numpy >= 1.19",
        "scipy >= 1.5.2",
        "matplotlib >= 3.2.2",
        "alive-progress >= 3.0.0"
    ],
    extras_require={
    },
)
