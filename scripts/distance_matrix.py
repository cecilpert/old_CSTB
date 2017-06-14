from ete3 import NCBITaxa
import numpy
import pickle

class Lineage: 
	def __init__(self):
		self.species="No specie"
		self.genus="No genus"	
		self.family="No family"
		self.order="No order"
		self.classe='No class'
		self.phylum="No phylum"

def dic_taxid(file): 
	f=open(file,'r')
	dic={}
	for l in f : 
		l_split=l.split('\t')
		ref=l_split[19].split('/')[-1].rstrip()
		taxid=l_split[5]
		dic[ref]=taxid

	f.close()	
	return dic

def create_lineage_objects(genome_file,dic_tax):
	ncbi=NCBITaxa()
	dic_lineage={}
	f=open(genome_file,'r')
	count=0
	for l in f: 
		lineage_object=Lineage()
		ref=l.split('\t')[1].rstrip().split('/')[-1]
		tax_ref=dic_tax[ref]
		lineage=ncbi.get_lineage(tax_ref)
		names=ncbi.get_taxid_translator(lineage)
		ranks=ncbi.get_rank(lineage)
		for i in ranks: 
			if ranks[i]=='species': 
				lineage_object.species=names[i]
			elif ranks[i]=='genus': 
				lineage_object.genus=names[i]
			elif ranks[i]=='family': 
				lineage_object.family=names[i]
			elif ranks[i]=='order': 
				lineage_object.order=names[i]
			elif ranks[i]=='class': 
				lineage_object.classe=names[i]	
			elif ranks[i]=='phylum': 
				lineage_object.phylum=names[i]					
		dic_lineage[ref]=(lineage_object,count)
		count+=1		
	f.close()
	return dic_lineage	

def create_matrix(dic_lineage): 
	matrix_tab=[]
	for a in dic_lineage: 
		sub_array=[]
		for b in dic_lineage: 
			if a==b: 
				sub_array.append(0)
			elif dic_lineage[a][0].species==dic_lineage[b][0].species: 
				sub_array.append(1)	
			elif dic_lineage[a][0].genus==dic_lineage[b][0].genus: 
				sub_array.append(2)	
			elif dic_lineage[a][0].family==dic_lineage[b][0].family: 
				sub_array.append(3)
			elif dic_lineage[a][0].order==dic_lineage[b][0].order: 
				sub_array.append(4)	
			elif dic_lineage[a][0].classe==dic_lineage[b][0].classe: 
				sub_array.append(5)		
			elif dic_lineage[a][0].phylum==dic_lineage[b][0].phylum: 
				sub_array.append(6)						
			else: 
				sub_array.append(10)	
		matrix_tab.append(sub_array)	
	matrix=numpy.matrix(matrix_tab)
	print(matrix)	

def distance_dic(dic_lineage): 
	dic={}
	for ref1 in dic_lineage: 
		dic[ref1]={}
		for ref2 in dic_lineage: 
			if not ref1==ref2: 
				if dic_lineage[ref1][0].species==dic_lineage[ref2][0].species: 
					dic[ref1][ref2]=1
				elif dic_lineage[ref1][0].genus==dic_lineage[ref2][0].genus: 
					dic[ref1][ref2]=2
				elif dic_lineage[ref1][0].family==dic_lineage[ref2][0].family: 
					dic[ref1][ref2]=3
				elif dic_lineage[ref1][0].order==dic_lineage[ref2][0].order: 
					dic[ref1][ref2]=4	
				elif dic_lineage[ref1][0].classe==dic_lineage[ref2][0].classe: 
					dic[ref1][ref2]=5		
				elif dic_lineage[ref1][0].phylum==dic_lineage[ref2][0].phylum: 
					dic[ref1][ref2]=6					
				else: 
					dic[ref1][ref2]=7
	return dic					


ref_file='more_genomes/assembly_summary_bacteria.txt'
genome_file='scripts/reference_genomes_ftp.txt'
dic_tax=dic_taxid(ref_file)
dic_lineage=create_lineage_objects(genome_file,dic_tax)
dist_dic=distance_dic(dic_lineage)
pickle.dump(dist_dic, open( "reference_genomes/distance_dic.pickle", "wb" ) )
