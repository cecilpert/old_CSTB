from flask import Flask, render_template, jsonify, request, send_file
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json

app=Flask(__name__)

PYTHON_INTERPRETER="python3"
ROOT_FOLDER="/Users/cecilehilpert/CSTB"
DATA_FOLDER="/Users/cecilehilpert/CSTB/reference_genomes"
CACHE_FOLDER="/Users/cecilehilpert/CSTB/tmp"

#PYTHON_INTERPRETER="python3"


# Configure we application port and adresss

def treat_arguments_allgenomes(): 
	#Genomes_IN and NOT_IN lists parsing.
	gi=request.args.get('gi',0).strip('[]').split(',')
	print('gi1',gi)
	gi=[i.strip('"\\n\\t') for i in gi]
	gi='"'+'&'.join(gi)+'"'	#Spaces escaped in argument using '"' character.
	gni=request.args.get('gni',0).strip('[]').split(',')
	gni=[i.strip('"\\n\\t') for i in gni]
	gni='"'+'&'.join(gni)+'"'
	print(gni)

	#Other parameters parsing

	pam=request.args.get('pam',0)
	pam=str(pam)
	pam=pam.replace('"','')

	sgrna_length=request.args.get('sgrna_length',0)
	sgrna_length=str(sgrna_length)
	sgrna_length=sgrna_length.replace('"','')

	command=PYTHON_INTERPRETER + " " + ROOT_FOLDER + "/scripts/allgenomes.py -cah " + CACHE_FOLDER + " -rfg " + DATA_FOLDER + " -gi " + gi + " -gni " + gni + " -pam " + pam + " -sl " + sgrna_length
	print(command)

	return command

def treat_results_all_genomes(command): 
	output=os.popen(command,'r')
	lines=output.readlines()
	for line in lines:
		if "Program terminated" in line:
			error_split=line.split('&')
			info=error_split[1].rstrip("\n]")
			return jsonify("Search yielded no results.",info)
	else:
		tag=lines[1].strip()
		with open(CACHE_FOLDER+'/'+tag+'/results.json','r') as f: 
			res=f.read()
		not_in=lines[0].strip()
		number_hits=lines[2].strip()
		return jsonify(res,not_in,tag,number_hits)		

def treat_arguments_specific_gene():
	seq=request.args.get('seq',0).strip('"')
	seq=str(seq)	## NB THIS REQUIRES PROPER FASTA FORMATTING IN JAVASCRIPT.

	#Gin and Gnotin treatment
	gin=request.args.get('gin',0).strip('[]').split(',')	##Relies on having no ',' in the organism names.
	gin=[genome.strip('"\\n\\t') for genome in gin]
	gin='"'+'&'.join(gin)+'"'
	gin=gin.rstrip('+"')+'"' #delete + at the end of the string if there is only one genome

	gnotin=request.args.get('gnotin',0).strip('[]').split(',')
	gnotin=[genome.strip('"\\n\\t') for genome in gnotin]
	gnotin='"'+'&'.join(gnotin)+'"'

	#Other parameters
	n=request.args.get('n',0)
	n=str(n)
	n=n.replace('"','')

	percent_id=request.args.get('pid',0)
	percent_id=str(percent_id)
	percent_id=percent_id.replace('"','')

	pam=request.args.get('pam',0)
	pam=str(pam)
	pam=pam.replace('"','')

	sgrna_length=request.args.get('sgrna_length',0)
	sgrna_length=str(sgrna_length)
	sgrna_length=sgrna_length.replace('"','')


	command=PYTHON_INTERPRETER + " " + ROOT_FOLDER + "/scripts/specificgene.py -cah " + CACHE_FOLDER + " -rfg " + DATA_FOLDER + " -seq " + seq + " -gi " + gin + " -gni " + gnotin + " -n " + n + " -ip " + percent_id + " -pam " + pam + " -sl " + sgrna_length 
	print(command)
	return command

def treat_results_specific_gene(command): 
	output=os.popen(command,'r')
	lines=output.readlines()	##List containing all print statements in specific gene script.
	#print(all_lines)
	for line in lines:
		if "Program terminated" in line:
			error_split=line.split('&')
			info=error_split[1].rstrip("\n]")
			return jsonify("Search yielded no results.",info)
	else:
		tag=lines[1].strip()
		with open(CACHE_FOLDER+'/'+tag+'/results.json','r') as f: 
			res=f.read()
		print(res)
		not_in=lines[0].strip()
		number_hits=lines[2].strip()
		number_on_gene=lines[3].strip()
		return jsonify(res,not_in,tag,number_hits,number_on_gene)

@app.route('/')
def rendu():	##HTML rendering upon url access
	return render_template('interface_layout.html')

@app.route('/allgenomes')
def ret():
	command=treat_arguments_allgenomes()
	result=treat_results_all_genomes(command)
	return result
	
@app.route('/specific_gene')
def ret_specific_gene():
	command=treat_arguments_specific_gene()
	result=treat_results_specific_gene(command)
	return result

@app.route('/download', defaults={'path': ''})
@app.route('/download/<path:path>')
def downloadResultFile(path):
	print ('Trying to server download request w/ key : %s' % path)
	filePath=CACHE_FOLDER+'/'+path+'/results_allgenome.txt'
	return send_file(filePath, mimetype='plaintext')	
	
