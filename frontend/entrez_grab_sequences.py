#!/usr/bin/python
import getopt
from Bio import Entrez,SeqIO
import sys

def get_entrez_ids(parameters):
	'''Get IDs for all species-hit-combinations out of database'''
	species_hash = {}
	species_errors = []
	
	# read species-names
	try:
		species_file = open(parameters["species_file"],"r")
		for line in species_file:
			species_hash[line.strip()] = 1
	except: 
		return [0,"error_in_species_file"]

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
					term=species+"[Orgn] AND "+parameters["search_term"],retmax=10**9)
				record = Entrez.read(handle)
				sys.stdout.write("\r")	
				hit_ids[species] = record["IdList"]
				handle.close()
				error_counter = 5
			except:
				error_counter += 1
		if error_counter == 3:
			print "\nThere was a problem to fetch the IDs for species "+species
			species_errors.append(species)

	print "\n3. Successfully got IDs for each species"
	return [hit_ids,species_errors]

def get_publications(parameters,hit_ids):
	'''get all publications out of pubmed-ids'''

	id_errors = []
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
				id_errors.append(single_id)
		print "\n5. Got all publications for all hits"
		return id_errors

def get_sequences(parameters,hit_ids):
	'''Get all sequences if database is nucleotide/protein'''
	out_file = open(parameters["out_file"],"w")
	sequences = []
	counter = 0
	id_errors = []
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
				id_errors.append(single_id)
	print "\n5. Successfully got all sequences for all hits"

	# write all sequences to outfile
	SeqIO.write(sequences,out_file,"fasta")
	out_file.close()	
	print "6. Wrote all sequences to ",parameters["out_file"]
	return id_errors			
