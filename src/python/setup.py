#from distutils.core import setup
import setuptools
from setuptools import setup
from os import path

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
  name = 'type3detect',         # How you named your package folder 
  packages = setuptools.find_packages(),#['lofarSun','lofarSun.IM','lofarSun.BF','lofarSun.BF.GUI'],   # Chose the same as "name"
  include_package_data=True,
  version = '0.0.3',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'tools to process the lofar solar data',   # Give a short description about your library
  author = 'Peijin',                   # Type in your name
  author_email = 'pjer1316@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/peijin94/type3detect',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/peijin94/type3detect/archive/refs/heads/master.zip',    
  keywords = ['LOFAR', 'Solar', 'radio'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'matplotlib',
          'sunpy',
          'astropy',
          'h5py',
          'scipy',
          'numpy'
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',     # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   # Again, pick a license,
    'Programming Language :: Python :: 3.9',
  ],
  long_description=long_description,
  long_description_content_type='text/markdown'
)
