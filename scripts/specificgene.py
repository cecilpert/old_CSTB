from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio import motifs
from Bio.Alphabet import IUPAC
from functionbase import reverse_complement,construct_dict_organism_assemblyref,build_expression,reverse_complement_mismatch,timestamp,cumulative_length,corrected_mismatch
import argparse, os, sys, pickle, time, re
import cProfile 


class Hit():
	def __init__(self,sequence,in_dic):
		self.sequence=sequence
		self.in_dic=in_dic
		self.notin_dic={}
		self.score=0

### METHODS FOR THE ARGUMENTS ###



def args_gestion():
	'''Take and treat arguments that user gives in command line'''
	parser = argparse.ArgumentParser(__file__, description="Specific gene program.")
	parser.add_argument("-seq", help="String query sequence", required=True)
	parser.add_argument("-gi", help="List with the name(s) of genome(s)", required=True)
	parser.add_argument("-n", help="The research will be done on the n first bases of the gene", required=True,type=int)
	parser.add_argument("-gni", help='List with the name(s) of not in genome(s)', required=True)
	parser.add_argument('-ip',
						help='identity percentage min for the research of homologous genes using blastn (default:70)', required=True,type=int)
	parser.add_argument('-mm',help='max mismatch between sgRNA', required=True,type=int)
	parser.add_argument('-pam',help='PAM motif for sgRNA', required=True)
	parser.add_argument('-sl',help='sgRNA length (without PAM motif)', required=True, type=int)
	parser.add_argument('-mmog',help='max mismatch between sgRNA and its origin genome', required=True, type=int)
	parser.add_argument('-mmnin',help='mismatch between sgRNA and excluded genome to be considered as absent', required=True, type=int)
	args = parser.parse_args()

	# Construct the list of origin and not in genome(s)
	dic_genome = construct_dict_organism_assemblyref()
	genome_list = args.gi.split('+')
	if args.gni:
		not_in_genome_list=args.gni.split('+')
	else:
		not_in_genome_list=[]

	return args.seq, genome_list, dic_genome, args.n, not_in_genome_list, args.ip, args.mm, args.pam, args.sl, args.mmog, args.mmnin



def eprint(*args, **kwargs):
	'''For printing only on the terminal'''
	print(*args, file=sys.stderr, **kwargs)


### METHODS FOR THE RESEARCH ###

def find_sgRNA_seq(seq,pam,size_sgRNA,organism):
	'''Find all sgRNA sequence using regular expression'''
	dict_seq={}
	pam=reverse_complement(pam)
	reg_exp = build_expression(pam) 
	indices = [m.start() for m in re.finditer('(?=' + reg_exp + ')', seq, re.I)]
	len_to_search=size_sgRNA+len(pam)
	for i in indices: 
		if i + len_to_search<=len(seq): 
			seq_sgRNA=seq[i:i+len_to_search]
			if seq_sgRNA not in dict_seq: 
				dict_seq[seq_sgRNA]=''
	return dict_seq


def blast_to_find_all_genes(query_seq, genomes, blast_path, path_reference_genomes, identity_percentage_min,dic_genome):
	'''Launch blastn between query sequence and reference genome and all in genomes (separately).
	For origin genome, will only save the genes position, not the subject sequence (the query sequence is conserved for origin genome)
	For other genomes IN, will take the subject sequence found by blast as sequence for the organism and the sgRNAs will be searched on it. It will also return gene's genomic position
	The hits must have a %id greater than the %id given by user 
	If there is several blast hits, we only take the first (as defined by Blast)
	'''

	list_genes=[]   

	size_gene=len(query_seq)
	file_query='./scripts/temp_blast.txt'
	with open(file_query,'w') as f:
		f.write(query_seq)
	f.close()    
	for i in range(len(genomes)): 
		eprint("Searching for sequence in",genomes[i])
		db_path = path_reference_genomes + "fasta/" + dic_genome[genomes[i]] + "_genomic.fna" #peut-etre mieux faire sur les cds parce que avec le génome complet pas d'infos sur nom gene, description etc...
		blast_command = blast_path + "blastn -db " + db_path + " -query " + file_query + " -outfmt 5"
		blast_output = os.popen(blast_command, 'r')
		blast_records = NCBIXML.parse(blast_output)
		for blast_record in blast_records:
			if blast_record.alignments == []:
				eprint("No blast hits for",genomes[i])
				print("Program terminated&No blast hits for "+genomes[i]+". Your sequence probably does not have homologs in this genome.")
				exit(1)
			else:
				first_alignment=blast_record.alignments[0]
				first_hsp=first_alignment.hsps[0]
				identity_percentage=first_hsp.identities/len(first_hsp.query)
				if identity_percentage < (identity_percentage_min/100): 
					eprint("Blast hit(s) do not meet stringency criterium: ID percentage. Program termination.")
					print("Program terminated& Blast hit in "+genomes[i]+" do not meet stringency criterium ID percentage. Maybe try with a lower ID percentage.")
					exit(1)
				else: 
					if i==0: #first genome, we want the query seq to search sgRNA 
						if first_hsp.sbjct_start > first_hsp.sbjct_end: #gene is on reverse strand, + treatment to have position of the query seq and not the position of the sbjct hit (both can be not exactly the same, and we want to search sgRNA in our query seq) 
							gene_start = first_hsp.sbjct_end - (size_gene - first_hsp.query_end) - 1
							gene_end = first_hsp.sbjct_start + first_hsp.query_start - 2
									
						else: #gene on forward strand, treatment to have position of the query seq 
							gene_start = first_hsp.sbjct_start - first_hsp.query_start
							gene_end = first_hsp.sbjct_end + (size_gene - first_hsp.query_end) - 1
								
						sublist=[query_seq,gene_start,gene_end]
						list_genes.append(sublist)

					else: #other in genomes, now we take the hit subject as seq for sgRNA searching  
						if first_hsp.sbjct_start > first_hsp.sbjct_end: #reverse strand 
							gene_start=first_hsp.sbjct_end-1 #blast is on 1-based and our program need 0-based 
							gene_end=first_hsp.sbjct_start-1
							strand='-'
						else:
							gene_start=first_hsp.sbjct_start-1
							gene_end=first_hsp.sbjct_end-1  
							strand='+'
						sublist=[first_hsp.sbjct,gene_start,gene_end]
						list_genes.append(sublist)   
	os.system('rm scripts/temp_blast.txt') #delete temporary file after use                                                   
	return list_genes   

										   

def search_in_origin_genome(dic_sgRNA_seq, genome_ref, bowtie_path, path_reference_genomes, gene_start,
							gene_end,pam,sgrna_length, mismatch_og):
	'''Align each sgRNA found in origin genome with this genome (using Bowtie) to found if it match at several places.
	Treat the results with treat_bowtie_results_in function'''
	out=open('tmp/sgRNA_in.fa','w')
	bowtie_index = path_reference_genomes + "index1/" + genome_ref
	for seq in dic_sgRNA_seq:
		out.write('>'+seq+'\n'+seq+'\n')
	out.close()    
	bowtie_command = bowtie_path + "bowtie --quiet -a -v "+str(mismatch_og)+" "+ bowtie_index + " -f tmp/sgRNA_in.fa tmp/bowtie_results.txt"  # align the sgRNA against its genome of origin
	#Bowtie options: -a=show all results, , -v=max number of tolerated mismatches to reference
	os.system(bowtie_command)
	dic_results=treat_bowtie_results_in(pam,sgrna_length,gene_start,gene_end)
	return dic_results
	

def treat_bowtie_results_in(pam,sgrna_length,gene_start,gene_end):
	'''Parse the result file from bowtie to create a dictionnary like key=sequence searched and value=coordinates found in the format (sequence_reference;strand;start;end;onGene or OffGene;mismatch)
	'''
	f=open('tmp/bowtie_results.txt','r')
	dic_results={}
	forbidden_pos=[]
	for i in range(len(pam)): 
		forbidden_pos.append(str(i))
	for l in f: 
		col=l.rstrip().split('\t')
		seq=col[0]
		strand=col[1]
		gene=col[2]
		start=int(col[3])
		end=start+len(pam)+int(sgrna_length)-1
		try: 
			mm=col[7]
		except:
			mm=False
		wrongmismatches=False
		if mm: 
			corrected_mm=''
			mm_list=mm.split(',')  
			for mm_el in mm_list: 
				if mm_el.split(':')[0] in forbidden_pos: 
					corrected_el='Mismatch in PAM'
				else:	
					corrected_el=corrected_mismatch(mm_el,sgrna_length+len(pam),strand)
				corrected_mm=corrected_el+','
			mm=corrected_mm.rstrip(',')	
		else: 
			mm='0'           
		if sgRNA_on_gene(gene_start,gene_end,start,end): 
			on_gene='OnGene'
		else: 
			on_gene='OffGene'   
		if seq not in dic_results:
			dic_results[seq]=['('+gene+';'+strand+';'+str(start)+';'+str(end)+';'+on_gene+';'+mm+')']
		else: 
			dic_results[seq].append('('+gene+';'+strand+';'+str(start)+';'+str(end)+';'+on_gene+';'+mm+')')           
	return dic_results

def treat_bowtie_results_not_in(pam,sgrna_length):
	'''Parse the result file from bowtie to create a dictionnary like key=sequence searched and value=coordinates found in the format (sequence_reference;strand;start;end;mismatch)
	'''
	f=open('tmp/bowtie_results_not_in.txt','r')
	dic_results={}
	forbidden_pos=[]
	for i in range(len(pam)): 
		forbidden_pos.append(str(i))
	for l in f: 
		col=l.rstrip().split('\t')
		seq=col[0]
		strand=col[1]
		gene=col[2]
		start=int(col[3])
		end=start+len(pam)+int(sgrna_length)-1
		try: 
			mm=col[7]
		except:
			mm=False
		wrongmismatches=False
		if mm: 
			corrected_mm=''
			mm_list=mm.split(',')  
			for mm_el in mm_list: 
				if mm_el.split(':')[0] in forbidden_pos: 
					corrected_el='Mismatch in PAM'
				else:	
					corrected_el=corrected_mismatch(mm_el,sgrna_length+len(pam),strand)
				corrected_mm=corrected_el+','
			mm=corrected_mm.rstrip(',')	
		else: 
			mm='0'              
		if seq not in dic_results:
			dic_results[seq]=['('+gene+';'+strand+';'+str(start)+';'+str(end)+';'+mm+')']
		else: 
			dic_results[seq].append('('+gene+';'+strand+';'+str(start)+';'+str(end)+';'+mm+')')           
	return dic_results


def sgRNA_on_gene(start_gene, end_gene, start_sgRNA, end_sgRNA):
	'''Check if a sgRNA is in a gene (with both genomics positions)'''
	on_gene = False
	if start_gene <= start_sgRNA and end_gene >= end_sgRNA:
		on_gene = True
	return on_gene


def find_common_sgRNA(genome_list,max_mismatches,dic_sgRNA,bowtie_path,path_reference_genomes,dic_genome,gene_list,pam,sgrna_length):
	'''Find the common sgRNAs between two or more sequences, with maximum mismatches given by user. '''
	for i in range(len(genome_list)):
		genome=genome_list[i]
		gene_start=gene_list[i][1]
		gene_end=gene_list[i][2] 
		new_dic={}
		bowtie_index = path_reference_genomes + "index1/" + dic_genome[genome]
		bowtie_command = bowtie_path + "bowtie --quiet -a -v "+str(max_mismatches)+" "+ bowtie_index + " -f tmp/sgRNA_in.fa tmp/bowtie_results.txt"  # align the sgRNA against its genome of origin
		os.system(bowtie_command)
		dic_results=treat_bowtie_results_in(pam,sgrna_length,gene_start,gene_end)
		out=open('tmp/sgRNA_in.fa','w')
		for seq in dic_results: 
			out.write('>'+seq+'\n'+seq+'\n')
			new_dic[seq]=dic_sgRNA[seq]
			new_dic[seq][genome]=dic_results[seq]
		out.close()    
		dic_sgRNA=new_dic
		if not dic_sgRNA:  # Empty list evaluates to False
			print("Program terminated& Adding organism "+genome+" gave an absence of common occurences. Maybe try to do the research with less organisms, with a further position in research or with a shorter sgRNA. \n")
			sys.exit(1)

	return(dic_sgRNA)    


def search_sgRNA_not_in(hitlist,not_in_genome_list,bowtie_path,mismatch,path_reference_genomes,dic_genome,pam,sgrna_length):
	'''For each excluded organism : 
	- do bowtie research with common sgRNA previously found
	- treat this results with treat_bowtie_results_not_in function, that give a result dictionnary like key=sequence and value=coordinates found by bowtie
	- complete the objects Hits with complete_hitlist function'''
	for genome in not_in_genome_list:
		bowtie_index=path_reference_genomes+'index1/'+dic_genome[genome]
		bowtie_command=bowtie_path + "bowtie --quiet -a -v "+str(mismatch)+" "+ bowtie_index + " -f tmp/sgRNA_in.fa tmp/bowtie_results_not_in.txt"  # align the sgRNA against its genome of origin
		os.system(bowtie_command)
		dic_results=treat_bowtie_results_not_in(pam,sgrna_length)
		hitlist=complete_hitlist(hitlist,dic_results,genome)
	return hitlist

def complete_hitlist(hitlist,dic_results,organism): 
	'''Complete all Hits objects with informations about not in research. Fill the attributåe notin_dic'''
	for hit in hitlist: 
		if hit.sequence in dic_results: 
			hit.notin_dic[organism]=dic_results[hit.sequence]
		else: 
			hit.notin_dic[organism]='Absent'	
	return hitlist			


def output_interface(hitlist,list_genome,list_notin_genome): 
	''''Construct json sorted output for research'''
	json_list=[]
	for hit in hitlist: 
		seq=reverse_complement(hit.sequence)
		ref_org=list_genome[0]
		ref_coordinates='%'.join(hit.in_dic[ref_org])
		if len(list_genome)>1: 
			otherorgs='%'.join(list_genome[1:])
			other_org_str='[{'
			for other_org in list_genome[1:]:
				coordinates='%'.join(hit.in_dic[other_org])
				other_org_str+='"'+other_org+'":"'+coordinates+'",'
			other_org_str=other_org_str.rstrip(',')+'}]'	
		else:
			otherorgs='No search'	
			other_org_str='""'
		
		if list_notin_genome==[]:
			not_in_str='"No search"'
			not_in_genomes_str=''	
		else: 
			not_in_str='[{'
			for org_notin in list_notin_genome: 
				if hit.notin_dic[org_notin]=='Absent':
					coordinates='Absent'
				else: 
					coordinates='%'.join(hit.notin_dic[org_notin])
				not_in_str+='"'+org_notin+'":"'+coordinates+'",'		 
			not_in_str=not_in_str.rstrip(',')+'}]'	
			not_in_genomes_str='%'.join(list_notin_genome)		

		jsonhit='{"refsequence":"'+seq+'","reforg":"'+ref_org+'","on_off":"'+ref_coordinates+'","otherorgs":"'+otherorgs+'","other_on_off":'+other_org_str+',"not_in":'+not_in_str+',"not_in_genomes":"'+not_in_genomes_str+'","score":"'+str(hit.score)+'"}'
	
		json_list.append(jsonhit)
	print("["+",".join([jstring for jstring in json_list])+"]")	


def write_result_file(hitlist,n,genome_list,identity_percentage_min,max_mismatch,sgrna_length,pam,not_in_genome_list):
	'''Write results in file'''
	new_tag=timestamp()
	print(new_tag)
	out=open('./scripts/static/results'+new_tag+'.txt','w')
	out.write('#genomic position format : (id chr/plasmid/complete genome:strand:start:end:OnGene or OffGene (target gene):mismatch with reference organism like given by Bowtie : 0 position is the first base).\n')     
	out.write("#mismatches : all mismatches are given compared to the reference seq pamNN..NNNNNN (not reverse NNNN...NNpam). Mismatches are given in 0-based (position 0 is the first base of the sequence) \n")
	out.write("@Parameters: ")
	if n==0: 
		out.write("search in all sequence, ")
	else:
		out.write("search for position 1 to "+str(n)+', ')
	out.write("identity percentage for homolog genes:"+str(identity_percentage_min)+", max mismatch:"+str(max_mismatch)+", sgRNA length:"+str(sgrna_length)+", PAM motif:"+pam+"\n")        
	out.write('@origin genome:'+genome_list[0]+' other targeted genomes:'+";".join(genome_list[1:]))  
	if not_in_genome_list:
		out.write(' excluded genomes:'+";".join(not_in_genome_list))
	out.write('\n')    
	out.write('sgRNA sequence\treference organism\tgenomic position(s)\t')
	for i in range(1,len(genome_list)):
		out.write("other organism\tgenomic position(s)\tmismatch(s)\t")
	for i in not_in_genome_list:    
		out.write("excluded organism\tpresence in excluded organism\t")
	out.write("score\n")

	for hit in hitlist: 
		reverse_sgrna=reverse_complement(hit.sequence)
		ref_org=genome_list[0]
		ref_coordinates='&'.join(hit.in_dic[ref_org])
		out.write(reverse_sgrna+'\t'+ref_org+'\t'+ref_coordinates)  
		if len(genome_list)>1:
			for other_org in genome_list[1:]: 
				coordinates='&'.join(hit.in_dic[other_org])
				out.write('\t'+other_org+'\t'+coordinates)  
		for org_notin in not_in_genome_list: 
			if hit.notin_dic[org_notin]!='Absent':
				coordinates='&'.join(hit.notin_dic[org_notin])
			else: 
				coordinates='Absent'
			out.write('\t'+org_notin+'\t'+coordinates)			
		out.write('\t'+str(hit.score))
		out.write('\n')    
	out.close()
             


def score_and_sort(hitlist):
	'''Score and sort the results of an analysis.
	For now, the score for 1 sgRNA is the sum of all mismatchs + 1 if it's present in an excluded organism (no matter if it's exact match or not)  
	Fill the attribute score of Hit objects.
	So best sgRNA is the one that have lowest score. 
	''' 
	for hit in hitlist:
		score=0
		for org in hit.in_dic: 
			for occurence in hit.in_dic[org]:
				mm=occurence.split(';')[-1].rstrip(')')
				if mm!='0': 
					score+=len(mm.split(','))
		for org_notin in hit.notin_dic: 
			if hit.notin_dic[org_notin]!='Absent':
				score+=1
		hit.score=score		
	sorted_hitlist=sorted(hitlist,key=lambda hit:hit.score)    
	return sorted_hitlist
 
			
def create_sgRNA_dic(organism,dic_results):
	'''Create a dictionary like key=sgRNA sequence and value=another dictionnary with key=organism and value=coordinates of the sequence in this organism'''
	dic_sgRNA={}
	for seq in dic_results: 
		dic_sgRNA[seq]={organism:dic_results[seq]}
	return dic_sgRNA     

def create_hitlist(dic_sgRNA): 
	'''Create hitlist that contains Hit objects. There's one Hit object for each sgRNA.
	Only the attributes sequences and in_dic are filled with this function'''
	hitlist=[]
	for seq in dic_sgRNA: 
		hit=Hit(seq,dic_sgRNA[seq])    
		hitlist.append(hit)
	return(hitlist)    

def do_search(query_seq, n, genome_list, dic_genome, not_in_genome_list,
			  identity_percentage_min,max_mismatch,pam,sgrna_length, mismatch_og, mismatch_notin):
	'''Launch the research with all parameters given by user. Principal function, it will call other functions to do the research'''

	##Results file 
	os.system('mkdir tmp')
	##Path for software used and database
	path_bowtie = "./bowtie-1.1.2/"
	path_reference_genomes = "./reference_genomes/"
	blast_path = "./ncbi-blast-2.5.0+/bin/"
	count = 0

	len_all_sgrna=sgrna_length+len(pam)
	all_sgRNA_list = []  ##This list will contain all_sgRNA objects which themselves contain occurence dict and sgRNA objects.

	##BLAST RESEARCH, to find all genes we need  
	list_genes=blast_to_find_all_genes(query_seq, genome_list, blast_path, path_reference_genomes, identity_percentage_min,dic_genome)

	first_gene=list_genes[0]
	first_genome=genome_list[0]
	gene_seq=str(first_gene[0])
	gene_start=first_gene[1]
	gene_end=first_gene[2]

	if n==0:
		dic_seq=find_sgRNA_seq(gene_seq,pam,sgrna_length,first_genome)
	else: 
		dic_seq=find_sgRNA_seq(gene_seq[:n], pam, sgrna_length,first_genome)    
	
	if dic_seq=={}: 
		eprint("No sgRNA's found in",genome_list[count],".Program termination.")
		print("Program terminated&No sgRNA's found in "+genome_list[count]+". Maybe try with a larger search position interval.")
		exit(1)   

	else: 
		dic_results=search_in_origin_genome(dic_seq,dic_genome[first_genome],path_bowtie,path_reference_genomes,gene_start,gene_end,pam,sgrna_length,mismatch_og)         
		dic_sgRNA=create_sgRNA_dic(first_genome,dic_results)
	eprint(str(len(dic_sgRNA))+' sgRNA in origin genome')    

	if len(genome_list) > 1: 
		eprint("Looking for common sgRNAs...")
		dic_sgRNA=find_common_sgRNA(genome_list[1:],max_mismatch,dic_sgRNA,path_bowtie,path_reference_genomes,dic_genome,list_genes[1:],pam,sgrna_length)
   
	eprint(str(len(dic_sgRNA))+' common sgRNA found in all genomes') 
	
	hitlist=create_hitlist(dic_sgRNA) 


	if not_in_genome_list:
		eprint("Working NOT IN...")
		hitlist=search_sgRNA_not_in(hitlist,not_in_genome_list,path_bowtie,mismatch_notin,path_reference_genomes,dic_genome,pam,sgrna_length)	

	sorted_hitlist=score_and_sort(hitlist)

	
	write_result_file(sorted_hitlist[:1000],n,genome_list,identity_percentage_min,max_mismatch,sgrna_length,pam,not_in_genome_list)

	output_interface(sorted_hitlist[:100],genome_list,not_in_genome_list)

	os.system('rm -r tmp')

def main():
	start_time=time.time()
	query_seq, genome_list, dic_genome, n, not_in_genome_list, identity_percentage_min, max_mismatch, pam, sgrna_length, mismatch_og, mismatch_notin  = args_gestion()
	do_search(query_seq, n, genome_list, dic_genome, not_in_genome_list,identity_percentage_min, max_mismatch, pam, sgrna_length, mismatch_og, mismatch_notin)
	end_time=time.time()
	total_time=end_time-start_time
	eprint('TIME',total_time)
	eprint('CUMULATIVE_LENGTH',cumulative_length(genome_list,dic_genome))

if __name__ == "__main__":
	main()
