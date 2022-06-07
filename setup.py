from pkg_resources import Requirement
from setuptools import setup
import os, sys
from dreem_nap import __version__

try:
    with open('requirements.txt') as f:
        requirements = f.read().splitlines()
except:
    with open('../requirements.txt') as f:
        requirements = f.read().splitlines()


setup(
   name='dreem_nap',
   version=__version__,
   description='New Analysis Pipeline: visualize the output of DREEM',
   author='Yves Martin des Taillades',
   author_email='yves@martin.yt',
   packages=['dreem_nap'],  #same as name
   package_dir={'dreem_nap': 'dreem_nap'},
   py_modules=[
       'dreem_nap/data_wrangler',
       'dreem_nap/database',
       'dreem_nap/plot',
       'dreem_nap/utils'
   ],
   include_package_data=True,
   install_requires=requirements, #external packages as dependencies
)