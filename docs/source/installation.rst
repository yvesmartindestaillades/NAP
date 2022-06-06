Installation Guide
------------------

Clone this repo on your computer:
::

   cd [path_to_repo]
   git clone https://github.com/yvesmartindestaillades/NAP

If you don't have Joe Yesselman's dreem installed on your computer, download a local copy.
::

   cd libs
   git clone https://github.com/jyesselm/dreem

Make sure that you have all of the Python libraries dependencies.
::

   sudo pip install -r requirements.txt

Validate the installation by running the test routine.
::

   python3 tests/DREEM_firebase_plots.py
