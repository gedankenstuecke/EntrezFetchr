from django.conf import settings
from celery.task import task
from entrez_grab_sequences import *
from helpers import *  
import os,time

@task()
def add(x, y):
	print x+y
	return x + y

@task()
def get_data(parameter_hash):
	Entrez.email = parameter_hash["email"]
	if check_binary(parameter_hash["species_file"]):
		send_file_error(parameter_hash)
		print "binary file"
		return "binary file"
	id_result = get_entrez_ids(parameter_hash)
	entrez_ids = id_result[0]
	species_errors = id_result[1]
	if species_errors == "error_in_species_file":
		send_file_error(parameter_hash)
		return "error with file"	
	if parameter_hash["database"] == "pubmed":
		id_errors = get_publications(parameter_hash,entrez_ids)
	else:
		id_errors = get_sequences(parameter_hash,entrez_ids)
	print "done"
	if species_errors != [] or id_errors != []:
		errors = True
	else:
		errors = False
	send_results(parameter_hash,errors,species_errors,id_errors)
	os.remove(parameter_hash["species_file"])
	delete_old_files(settings.BASE_DIRECTORY+'static/downloads/',settings.TIME_THRESHOLD)
	print "removed input-file"
	return "done"

def delete_old_files(directory,threshold):
	now = time.time()
	threshold_time = now - 60*60*24*threshold
	for single_file in os.listdir(directory):
		if single_file != "readme.txt":
			creation_time = os.path.getctime(directory+single_file)
			if creation_time < threshold_time:
				print "deleted file"
				os.remove(directory+single_file)
			else:
				print "file is not old enough"
	return "done"
