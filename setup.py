from pkg_resources import Requirement
from setuptools import setup

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
   name='NAP',
   version='1.0',
   description='New Analysis Pipeline: visualize the output of DREEM',
   author='Yves Martin des Taillades',
   author_email='yves@martin.yt',
   packages=['NAP'],  #same as name
   install_requires=requirements, #external packages as dependencies
)