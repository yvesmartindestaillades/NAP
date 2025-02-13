from pkg_resources import Requirement
from setuptools import setup
import os, sys
from dreem_nap import __version__
import dreem_nap
import sys

try:
    with open('requirements.txt') as f:
        requirements = f.read().splitlines()
except:
    with open('../requirements.txt') as f:
        requirements = f.read().splitlines()

PYTHON_VERSION = (3,10)

if sys.version_info < PYTHON_VERSION:
    sys.exit(f"Python >= {PYTHON_VERSION[0]}.{PYTHON_VERSION[1]} required.")


setup(
   name='dreem_nap',
   version=__version__,
   license="MIT",
   description='New Analysis Pipeline: visualize the output of DREEM',
   author='Yves Martin des Taillades',
   author_email='yves@martin.yt',
   long_description= 'TODO',
   packages=['dreem_nap'],  #same as name
   package_dir={'dreem_nap': 'dreem_nap'},
   py_modules=[
       'dreem_nap/data_wrangler',
       'dreem_nap/database',
       'dreem_nap/plot',
       'dreem_nap/util', 
       'dreem_nap/data_manip',
       'dreem_nap/wrapper',
       'dreem_nap/study',
       'dreem_nap/data',
       'dreem_nap/clustering',
       'dreem_nap/deltaG',
   ],
   include_package_data=True,
   install_requires=requirements, #external packages as dependencies
   
)