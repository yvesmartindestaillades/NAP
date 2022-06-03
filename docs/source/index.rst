.. NAP documentation master file, created by
   sphinx-quickstart on Thu Jun  2 17:30:09 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to NAP's documentation!
===============================

New Analysis Pipeline (NAP) is written to analyse DREEM output.

Installation Guide
------------------

Clone this repo on your computer:
::

   cd [path_to_repo]
   git clone https://github.com/yvesmartindestaillades/NAP

If you don't have Joe Yesselman's dreem installed on your computer, download a local copy:
::

   cd libs
   git clone https://github.com/jyesselm/dreem

Make sure that you have all of the Python libraries dependencies
::

   sudo pip3 install pandas pickle-mixin firebase_admin numpy matplotlib python-string-utils datetime seaborn

.. toctree::
   usage

.. note::

   This project is under active development.



NAP has its documentation hosted on Read the Docs.