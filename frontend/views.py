from django.conf import settings
from django.shortcuts import render_to_response
from django.http import HttpResponse
from django.template import RequestContext
from frontend.forms import *
from datetime import datetime
from tasks import *

# Create your views here.

def index(request):
	if request.method == "POST":
		form = SubmitForm(request.POST,request.FILES)
		if form.is_valid():
			parameter_hash = {}
			parameter_hash["email"] = form.cleaned_data["email"]
			parameter_hash["database"] = form.cleaned_data["database"]
			parameter_hash["search_term"] = form.cleaned_data["search_term"]
			filename = str(datetime.now()).replace(":","-").replace(".","-").replace(" ","-")
			file_suffix = form.cleaned_data["database"]
			filename = settings.BASE_DIRECTORY+'static/uploads/' + filename + "_" + file_suffix + ".txt"
			with open(filename,"wb+") as destination:
				for chunk in request.FILES["species_file"]:
					destination.write(chunk)
			parameter_hash["species_file"] = filename
			if parameter_hash["database"] == "pubmed":
				parameter_hash["out_file"] = filename.replace(".txt",".csv").replace("static/uploads","static/downloads")
			else:
				parameter_hash["out_file"] = filename.replace(".txt",".fasta").replace("static/uploads","static/downloads")
			get_data.delay(parameter_hash)
			return render_to_response("frontend/send.html",
					{"database":parameter_hash["database"],
					 "email": parameter_hash["email"],
					 "search_term": parameter_hash["search_term"]},
					context_instance=RequestContext(request))
		else:
			return render_to_response("frontend/index.html",{"form":form},context_instance=RequestContext(request))
	else:
		form = SubmitForm()
		return render_to_response("frontend/index.html",{"form":form},context_instance=RequestContext(request))		
