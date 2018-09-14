
import setuptools
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='rastersmith',
      version='0.0.2',
      description='RasterSmith is a package to preprocess different NASA Earth observing satellite data products into common resolution, spatial reference, and format for easy analysis and processing across sensors.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='http://github.com/kmarkert/rastersmith',
      packages=setuptools.find_packages(),
      author='Kel Markert',
      author_email='kel.markert@gmail.com',
      license='MIT',
      zip_safe=False)
