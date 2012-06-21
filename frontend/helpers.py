from django.conf import settings
from django.core.mail import EmailMessage
from django.template.loader import get_template
from django.template import RequestContext,Context

def send_results(parameter_hash,errors,species_errors,id_errors):
	plaintext = get_template("email.txt")
	context = Context({"download":parameter_hash["out_file"].replace(settings.BASE_DIRECTORY+"static",""),
				"errors":errors,
				"species_errors":species_errors,
				"id_errors":id_errors,
				"database":parameter_hash["database"],
				"search_term":parameter_hash["search_term"]})

	plain_content = plaintext.render(context)
	email = EmailMessage("Your Entrez-Results are ready",
				plain_content,
				settings.EMAIL_HOST_USER,
				[parameter_hash["email"]])
	email.send()
	print "send email to "+parameter_hash["email"]

def send_file_error(parameter_hash):
	plaintext = get_template("email_file_error.txt")
	context = Context({"database":parameter_hash["database"],
			"search_term":parameter_hash["search_term"]})
	plain_content = plaintext.render(context)
	email = EmailMessage("EntrezFetch: A problem with your data occured", 
				plain_content, 
				settings.EMAIL_HOST_USER,[parameter_hash["email"]])
	email.send()
	print "send error-message to "+parameter_hash["email"]

def send_length_error(parameter_hash):
	plaintext = get_template("email_length_error.txt")
	context = Context({"database":parameter_hash["database"],
			"search_term":parameter_hash["search_term"],
			"threshold":settings.RESULTS_THRESHOLD})
	plain_content = plaintext.render(context)
	email = EmailMessage("EntrezFetch: Too many results", 
				plain_content, 
				settings.EMAIL_HOST_USER,[parameter_hash["email"]])
	email.send()
	print "send too-many-results-error to "+parameter_hash["email"]

def check_binary(file_name):
	file_handle = open(file_name, 'rb')
	try:
        	CHUNKSIZE = 1024
		while 1:
			chunk = file_handle.read(CHUNKSIZE)
			if '\0' in chunk:
				return True
			if len(chunk) < CHUNKSIZE:
				break 
	finally:
		file_handle.close()
	return False

