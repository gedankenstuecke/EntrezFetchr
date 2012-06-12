#!/usr/bin/python
import getopt
from Bio import Entrez,SeqIO
import sys

def get_arguments(arguments):
	'''Check if all parameters are present'''
	out_args = {}

	if len(arguments) == 1:
		print "use -h to get help on the available options and usage"
		sys.exit()
	else:
		arguments.pop(0)
		optlist, arguments = getopt.getopt(arguments, 'hi:s:d:o:e:')
		for arg,opt in optlist:
			if arg == "-i":
				out_args["species_file"] = opt
			elif arg == "-s":
				out_args["search_term"] = opt
			elif arg == "-d":
				out_args["database"] = opt
			elif arg == "-o":
				out_args["out_file"] = opt
			elif arg == "-e":
				out_args["email"] = opt
			elif arg == "-h":
				print "Welcome to the Entrez-Fetch-Script"
				print "---"
				print "This script allows you to search any Entrez-database for"
				print "sequences and publications of a multitude of species and a given search term"
				print "at the same time"
				print "You need to specify the following options:"
				print "-i: A file of species you are interested in, one species per line"
				print "-s: The search term you are interested in"
				print "-d: The database you want to search ('nucleotide' for xNA-sequences, 'protein' for AA-sequences, 'pubmed' for publications"
				print "-o: The outfile where sequences will be saved"
				print "-e: Your eMail-address (NCBI requires this)"
				print "-h (optional): displays this helpful help-message"
				sys.exit()
		if len(out_args) != 5:
			print "You need to specify all options!"
			print "Use -h to get help"
			sys.exit()
		if out_args["database"] != "pubmed" and out_args["database"] != "nucleotide" and out_args["database"] != "protein":
			print "Your database-choice is not supported"
			print "For now 'pubmed','nucleotide' and 'protein' are supported."
			sys.exit()
		print "1. Successfully read parameters"
		return out_args

def get_entrez_ids(parameters):
	'''Get IDs for all species-hit-combinations out of database'''
	species_hash = {}
	
	# read species-names
	species_file = open(parameters["species_file"],"r")
	for line in species_file:
		species_hash[line.strip()] = 1

	hit_ids = {}
	counter = 0.0
	print "2. Getting IDs for each species in species file:"
	
	# iterate over all species and grab database-IDs of all hits
	for species in species_hash:
		counter += 1
		progress = int(50.0*counter/len(species_hash))* "#"
		progress += (50 - len(progress))* " "
		sys.stderr.write("\r" + "0%"+ progress+"100%")
		error_counter = 0
		while error_counter < 3:
			try:
				handle = Entrez.esearch(db=parameters["database"],
					term=species+"[Orgn] AND "+parameters["search_term"])
				record = Entrez.read(handle)
				sys.stdout.write("\r")	
				hit_ids[species] = record["IdList"]
				handle.close()
				error_counter = 5
			except:
				error_counter += 1
		if error_counter == 3:
			print "\nThere was a problem to fetch the IDs for species "+species

	print "\n3. Successfully got IDs for each species"
	return hit_ids

def get_publications(parameters,hit_ids):
	'''get all publications out of pubmed-ids'''

	out_file = open(parameters["out_file"],"w")
	out_file.write("PMID\tAuthor\tTitle\tJournal\n")
	counter = 0
	print "4. Getting publication-details for all hits"
	for species,ids in hit_ids.items():
		counter += 1
		progress = int(50.0*counter/len(hit_ids))*"#"
		progress += (50-len(progress))* " "		
		sys.stderr.write("\r"+"0%"+progress+"100%")
		for single_id in ids:
			error_counter = 0
			while error_counter < 3:
				try:
					handle = Entrez.efetch(db=parameters["database"],
								id=single_id,
								retmode="xml")
					result = Entrez.read(handle)
					out_line = single_id+"\t"
					out_line += result[0]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"]+"\t"
					out_line += result[0]["MedlineCitation"]["Article"]["ArticleTitle"]+"\t"
					out_line += result[0]["MedlineCitation"]["Article"]["Journal"]["Title"]+"\n"
					out_file.write(out_line.encode('UTF-8'))
					error_counter = 5
					handle.close()
				except:
					error_counter += 1
			if error_counter == 3:
				print "\nThere was a problem to fetch the publication "+str(single_id)
		print "\n5. Got all publications for all hits"

def get_sequences(parameters,hit_ids):
	'''Get all sequences if database is nucleotide/protein'''
	out_file = open(parameters["out_file"],"w")
	sequences = []
	counter = 0
	print "4. Getting sequences for all hits"	

	# iterate over all hit-ids and get sequences
	for species,ids in hit_ids.items():
		counter += 1
		progress = int(50.0*counter/len(hit_ids))* "#"
		progress += (50- len(progress))* " "
		sys.stderr.write("\r" + "0%" + progress+"100%") 
		for single_id in ids:
			error_counter = 0
			while error_counter < 3:
				try:
					handle = Entrez.efetch(db=parameters["database"],
						id=single_id,
						rettype="gb",
						retmode="text")
					sequence = SeqIO.read(handle,"genbank")
					sequences.append(sequence)
					handle.close()
					error_counter = 4
				except:
					error_counter += 1
			if error_counter == 3:
				print "\nThere was an Error fetching the ID "+str(single_id)
	print "\n5. Successfully got all sequences for all hits"

	# write all sequences to outfile
	SeqIO.write(sequences,out_file,"fasta")
	out_file.close()	
	print "6. Wrote all sequences to ",parameters["out_file"]			

def main():
	print ""
	print "#################################################################"
	print "#  Will now try to get all hits matching species & search term  #"
	print "#################################################################"
	print ""

	parameters = get_arguments(sys.argv)
	Entrez.email = parameters["email"]
	entrez_ids = get_entrez_ids(parameters)
	if parameters["database"] == "pubmed":
		get_publications(parameters,entrez_ids)
	elif parameters["database"] == "nucleotide" or parameters["database"] == "protein":
		get_sequences(parameters,entrez_ids)
	else:
		print "4. Unkown database. Choices for -d are: pubmed, nucleotide, protein"
main()
