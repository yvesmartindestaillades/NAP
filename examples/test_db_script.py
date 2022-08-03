
from dreem_nap import database
from firebase_admin import credentials
import firebase_admin
from firebase_admin import db
import json

firebase_credentials_file = '../data/credentials_firebase.json'
with open(firebase_credentials_file) as file:
	firebase_credentials = json.load(file)

database.connect(firebase_credentials)
database.connect(firebase_credentials)

default_app = firebase_admin.initialize_app(cred, {
    'databaseURL': f"https://{cred['project_id']}-default-rtdb.firebaseio.com/"
    }) 

print(f"https://{cred['project_id']}-default-rtdb.firebaseio.com/")
