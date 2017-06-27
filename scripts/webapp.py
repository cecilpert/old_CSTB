from flask import Flask, render_template, jsonify, request, send_file
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json

app=Flask(__name__)

PYTHON_INTERPRETER="python3"
ROOT_FOLDER="/Users/cecilehilpert/CSTB"
DATA_FOLDER="/Users/cecilehilpert/CSTB"
CACHE_FOLDER="/Users/cecilehilpert/CSTB/tmp"

#PYTHON_INTERPRETER="python3"


# Configure we application port and adresss



@app.route('/')
def rendu():	##HTML rendering upon url access
	return render_template('interface_layout.html')

@app.route('/allgenomes')
def ret():

	#Genomes_IN and NOT_IN lists parsing.
	gi=request.args.get('gi',0).strip('[]').split(',')
	print('gi1',gi)
	gi=[i.strip('"\\n\\t') for i in gi]
	gi='"'+'+'.join(gi)+'"'	#Spaces escaped in argument using '"' character.
	gni=request.args.get('gni',0).strip('[]').split(',')
	gni=[i.strip('"\\n\\t') for i in gni]
	gni='"'+'+'.join(gni)+'"'
	print(gni)

	#Other parameters parsing

	pam=request.args.get('pam',0)
	pam=str(pam)
	pam=pam.replace('"','')

	sgrna_length=request.args.get('sgrna_length',0)
	sgrna_length=str(sgrna_length)
	sgrna_length=sgrna_length.replace('"','')

	command=PYTHON_INTERPRETER + " " + ROOT_FOLDER + "/scripts/allgenomes.py -cah " + CACHE_FOLDER + " -rfg " + DATA_FOLDER + "/reference_genomes -gi " + gi + " -gni " + gni + " -pam " + pam + " -sl " + sgrna_length
	print(command)
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
		print(res)
		not_in=lines[0].strip()
		number_hits=lines[2].strip()
		return jsonify(res,not_in,tag,number_hits)


@app.route('/download', defaults={'path': ''})
@app.route('/download/<path:path>')
def downloadResultFile(path):
	print ('Trying to server download request w/ key : %s' % path)
	filePath=CACHE_FOLDER+'/'+path+'/results_allgenome.txt'
	return send_file(filePath, mimetype='plaintext')


@app.route('/specific_gene')
def ret_specific_gene():

	#Seq treatment
	seq=request.args.get('seq',0).strip('"')
	seq=str(seq)	## NB THIS REQUIRES PROPER FASTA FORMATTING IN JAVASCRIPT.

	#Gin and Gnotin treatment
	gref=request.args.get('gref',0).strip('[]').replace('"','')
	gin=request.args.get('gin',0).strip('[]').split(',')	##Relies on having no ',' in the organism names.
	gin=[genome.strip('"\\n\\t') for genome in gin]
	gin.insert(0,gref)##Places genome of origin at first position of genomes in list.
	gin='"'+'+'.join(gin)+'"'
	gin=gin.rstrip('+"')+'"' #delete + at the end of the string if there is only one genome


	gnotin=request.args.get('gnotin',0).strip('[]').split(',')
	gnotin=[genome.strip('"\\n\\t') for genome in gnotin]
	gnotin='"'+'+'.join(gnotin)+'"'


	#Other parameters
	n=request.args.get('n',0)
	n=str(n)
	n=n.replace('"','')

	percent_id=request.args.get('pid',0)
	percent_id=str(percent_id)
	percent_id=percent_id.replace('"','')

	max_mismatch=request.args.get('max_mismatch',0)
	max_mismatch=str(max_mismatch)
	max_mismatch=max_mismatch.replace('"','')

	pam=request.args.get('pam',0)
	pam=str(pam)
	pam=pam.replace('"','')

	sgrna_length=request.args.get('sgrna_length',0)
	sgrna_length=str(sgrna_length)
	sgrna_length=sgrna_length.replace('"','')

	mm_og=request.args.get('max_mm_og',0)
	mm_og=str(mm_og)
	mm_og=mm_og.replace('"','')

	mm_nin=request.args.get('max_mm_notin',0)
	mm_nin=str(mm_nin)
	mm_nin=mm_nin.replace('"','')

	command=PYTHON_INTERPRETER + " " + ROOT_FOLDER + "/scripts/specificgene.py -seq " + seq + " -gi " + gin + " -gni " + gnotin + " -n " + n + " -ip " + percent_id + " -mm " + max_mismatch + " -pam " + pam + " -sl " + sgrna_length + " -mmog " + mm_og + " -mmnin " + mm_nin
	print(command)
	output=os.popen(command,'r')
	all_lines=output.readlines()	##List containing all print statements in specific gene script.
	#print(all_lines)
	for line in all_lines:
		if "Program terminated" in line:
			error_split=line.split('&')
			info=error_split[1].rstrip("\n]")
			return jsonify("Search yielded no results.",info)
	else:
		tag=all_lines[-2].rstrip()
		json_string=all_lines[-1].rstrip('\n')
		print(json_string)
		return jsonify(json_string,tag,'')
