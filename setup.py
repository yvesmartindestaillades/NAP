from pkg_resources import Requirement
from setuptools import setup
import os

try:
    with open('requirements.txt') as f:
        requirements = f.read().splitlines()
except:
    with open('../requirements.txt') as f:
        requirements = f.read().splitlines()
version = '1.0'

setup(
   name='dreem_nap',
   version=version,
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