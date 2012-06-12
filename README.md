# EntrezFetchr


A small Django project which uses BioPython to fetch Sequences &amp; Publications out of NCBIs Entrez databases.
Users can upload a list of species they are interested in and provide additional search terms to get back FASTA-formatted sequences out of the nucleotide/protein-database or get a csv-list of publications matching their organisms of interest and the search term. 

## Dependencies
The application makes use of and requires:
* Django (https://www.djangoproject.com/)
* Celery (http://celeryproject.org/)
* RabbitMQ (http://www.rabbitmq.com/)
* Biopython (http://biopython.org/)

## How-To
* Adjust the settings.py to your needs (location of the project, credentials for your SMTP-server)
* Run "./manage.py syncdb" to create the SQLite-database for storing celery-results.
* Start your RabbitMQ-server
* Start a Celery-worker: "./manage.py celeryd"
* Fire up the development-server "./manage.py runserver 0.0.0.0:8000"
* Point your browser to localhost:8000 and view the result

## Command-Line-Version
If you do not need the web-stuff and just want to run your search queries from the shell you can use EntrezFetchr.py which is located in /scripts. To use it just pass the script the parameters:
* -i for the species-input-file 
* -s for the search term you are interested in
* -d for the database you want to crawl, right now limited to pubmed, nucleotide and protein
* -e for your email-address, as Entrez requires you to provide contact details to use their database
* -o for the output file

Reminder: You still need Python & Biopython to run the script, but you can skip the Django/RabbitMQ/Celery-dependencies

## Design-Stuff
The CSS is Bootstrap, from Twitter (http://twitter.github.com/bootstrap/)
The icons are Glyphicons Free (http://glyphicons.com/) 

## Contact
If you encounter any bugs or have feature requests feel free to contact me at bgreshake@googlemail.com or use the issues at GitHub: https://github.com/gedankenstuecke/EntrezFetchr
