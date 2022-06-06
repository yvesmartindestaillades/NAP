from pkg_resources import Requirement
from setuptools import setup
import os

with open('../requirements.txt') as f:
    requirements = f.read().splitlines()

version = '1.0'

setup(
   name='NAP',
   version=version,
   description='New Analysis Pipeline: visualize the output of DREEM',
   author='Yves Martin des Taillades',
   author_email='yves@martin.yt',
   packages=['NAP'],  #same as name
   package_dir={'NAP': 'NAP'},
   py_modules=[
       'NAP/data_wrangler',
       'NAP/database',
       'NAP/plot',
       'NAP/utils'
   ],
   include_package_data=True,
   install_requires=requirements, #external packages as dependencies
)