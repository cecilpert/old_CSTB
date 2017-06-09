from Bio.Seq import Seq
from Bio import SeqIO
from functionbase import *
import time,argparse,os,sys,re,random
import cProfile


"""
---AllGenomes option: search for sgRNA constructs in whole genomes---

*Object breakdown:
-sgRNA: holds coordinates and strand it is found at.
-hit: an sgRNA that is in common to the genomes specified in genomes_IN. It contains the sequence, a genomes_Dict of the form {organism:[occurences]},
a NOTIN_Dict of the form: {organism: [mismatches/exactmatches]/"nomatches" ('/' denote OR)}, and a score in tuple format (IN score,NOTIN score)


"""


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)     ##Found at http://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python


class Hit():
    def __init__(self,sequence,genomes_Dict):
        self.sequence=sequence
        self.genomes_Dict=genomes_Dict
        self.score=0   ##Holds score as tuple of occurence score (on genomes_IN) and NOTIN score 


def args_gestion(dict_organism_code):
    '''Take and treat arguments that user gives in command line'''
    ##Argparsing
    parser=argparse.ArgumentParser(description="Allgenomes program.")
    parser.add_argument("-gi",metavar="<str>",help="The organisms to search inclusion in.",required=True)
    parser.add_argument("-gni",metavar="<str>",help="The organisms to search exclusion from",required=True)
    parser.add_argument("-sl",metavar="<int>",help="The length of the sgRNA, excluding PAM")
    parser.add_argument("-pam",metavar="<str>",help="The PAM motif",required=True)
    args=parser.parse_args()
    ##
    eprint('IN',args.gi)
    eprint('NOT IN',args.gni)
    organisms_selected=args.gi.split('+')
    organisms_excluded=args.gni.split('+')
    organisms_selected=[i for i in organisms_selected if i in dict_organism_code]
    organisms_excluded=[i for i in organisms_excluded if i in dict_organism_code]
    PAM=reverse_complement(args.pam)    ##The PAM is given as 'NGG'; search is as 'CCN'
    non_PAM_motif_length=int(args.sl)

    return organisms_selected,organisms_excluded,PAM,non_PAM_motif_length

def find_sgRNA_seq(seq,pam):
    '''
    Uses Regular expression matching of the PAM motif to the reference genome to get the start 
    positions (0-based) of each match
    '''
    list_seq=[]
    reg_exp = build_expression(pam) 
    indices = [m.start() for m in re.finditer('(?=' + reg_exp + ')', seq, re.I)]
    return indices
    

def construct_in(fasta_path,organism,organism_code,PAM,non_PAM_motif_length):
    '''
    Same function as construct_seq_dict except the dictionnary will not be like value=list of coordinates. It add the information about organism with dictionnary values like : dictionnary with key=organism and value=list of coordinates in this organism
    '''
    start=time.time()
    fasta_file=fasta_path +'/' + organism_code +'_genomic.fna'
    genome_seqrecord=next(SeqIO.parse(fasta_file, 'fasta'))

    genome_length=len(genome_seqrecord)
    eprint('length',genome_length)

    genome_seq=str(genome_seqrecord.seq)
    
    eprint(PAM)
    eprint(reverse_complement(PAM))
    seq_list_forward=find_sgRNA_seq(genome_seq,PAM)
    seq_list_reverse=find_sgRNA_seq(genome_seq,reverse_complement(PAM))

    seq_dict={}

    for indice in seq_list_forward: 
        end=indice+len(PAM)+non_PAM_motif_length
        if not end > genome_length:
            seq=str(genome_seqrecord.seq[indice:end].reverse_complement())
            if seq not in seq_dict:
                seq_dict[seq]={organism:[]}
            seq_dict[seq][organism].append('+('+str(indice+1)+','+str(end)+')')
            
    for indice in seq_list_reverse: 
        end=indice+len(PAM)
        start=indice-non_PAM_motif_length
        if not start<=0:
            seq=genome_seq[start:end]
            if seq not in seq_dict:
                seq_dict[seq]={organism:[]}   
            seq_dict[seq][organism].append('-('+str(start+1)+','+str(end)+')')
    end=time.time()
    return(seq_dict) 


def sort_genomes(list_genomes,fasta_path,dict_org_code): 
    '''Sort genomes by ascending size'''
    tmp_list=[]
    for genome in list_genomes: 
        fasta_genome=next(SeqIO.parse(fasta_path +'/' + dict_org_code[genome] +'_genomic.fna', 'fasta'))
        tmp_list.append((len(fasta_genome.seq),genome))
    genomes_sorted=[i[1] for i in sorted(tmp_list,key=lambda genome:genome[0])]   ##Sort by ascending size
    return(genomes_sorted)

def sort_genomes_desc(list_genomes,fasta_path,dict_org_code): 
    '''Sort genomes by ascending size'''
    tmp_list=[]
    for genome in list_genomes: 
        fasta_genome=next(SeqIO.parse(fasta_path +'/' + dict_org_code[genome] +'_genomic.fna', 'fasta'))
        tmp_list.append((len(fasta_genome.seq),genome))
    genomes_sorted=[i[1] for i in sorted(tmp_list,key=lambda genome:genome[0],reverse=True)]   ##Sort by ascending size
    return(genomes_sorted)    

def construct_hitlist(dict_seq):
    '''
    Will construct an object Hit for each sequence in dictionnary, and store all Hits in a list.
    These function only fill the attributes sequence and genomes_Dict of the object 
    '''
    hits_list=[]
    count=0
    for seq in dict_seq:
        count+=1
        new_hit=Hit(seq,dict_seq[seq])
        hits_list.append(new_hit)  
    return(hits_list)


def write_to_fasta(dic_seq): 
    '''Write sequences in dic_seq in a fasta file'''
    out=open('tmp/sgRNA.fa','w')
    count=0
    for seq in dic_seq: 
        count+=1
        out.write('>'+seq+'\n'+seq+'\n')
    out.close()            

def add_notin(organism_code,indexs_path,dic_seq): 
    bowtie_command='bowtie2 -x'+indexs_path+'2/'+organism_code+' -f tmp/sgRNA.fa -S tmp/results_bowtie.sam -L 13 -a --quiet'
    os.system(bowtie_command)
    new_dic=treat_bowtie_not_in(dic_seq)
    return new_dic

def add_in(bowtie_path,organism_code,indexs_path,dic_seq,genome,len_sgrna):  
    bowtie_command='bowtie2 -x '+indexs_path+'2/'+organism_code+' -f tmp/sgRNA.fa -S tmp/results_bowtie.sam -L 13 -a --quiet'  
    os.system(bowtie_command)
    new_dic=treat_bowtie_in(dic_seq,genome,len_sgrna)   
    return new_dic

def treat_bowtie_not_in(dic_seq): 
    res=open('tmp/results_bowtie.sam','r')
    count=0
    new_dic={}
    for l in res: 
        if l[0]!='@': 
            if l.split('\t')[2]=='*':
                seq=l.split('\t')[0]
                new_dic[seq]=dic_seq[seq]
    return new_dic  

def treat_bowtie_in(dic_seq,genome,len_sgrna): 
    res=open('tmp/results_bowtie.sam','r')
    count=0
    new_dic={}
    for l in res: 
        if l[0]!='@': 
            if l.split('\t')[2]!='*':
                l_split=l.split('\t')
                cigar=l_split[5]
                if cigar=='23M': 
                    mm=l_split[-2]
                    if mm.split(':')[-1]=='23': 
                        seq=l_split[0]
                        if seq not in new_dic: 
                            new_dic[seq]=dic_seq[seq]
                        if genome not in new_dic[seq]: 
                            new_dic[seq][genome]=[]    
                        seq_align=l_split[9]
                        if seq != seq_align: 
                            strand='+'
                        else:
                            strand='-'  
                        start=l_split[3]
                        end=int(start)+len_sgrna-1
                        coord=strand+'('+start+','+str(end)+')'
                        new_dic[seq][genome].append(coord)
    return new_dic                    
   
def Scorage_triage(hitlist):
    '''
    Scoring of the hits found, where positive scores mean stronger sgRNA constructs.
    Complete Hit objects with score and sort the hits with this scores. 
    '''
    for hit in hitlist: 
        for genome in hit.genomes_Dict: 
            hit_score+=len(hit.genomes_Dict[genome])
    sorted_hitlist=sorted(hitlist,key=lambda hit:hit.score,reverse=True)
    return(sorted_hitlist)
  
             
def write_to_file(genomes_IN,genomes_NOT_IN,hit_list,PAM,non_PAM_motif_length): 
    '''Write results in a file. 
    The file is a tabulated file, with first column=sgRNA sequence, then one column for each included organisms with list of positions of the sequence in each, and finally one column for each excluded organism with informations about alignment of the sequence with this genomes organisms.
    '''
    new_tag=timestamp()
    print(new_tag)
    output=open('./scripts/static/results'+new_tag+'.txt','w')
    not_in=True
    gi=','.join(genomes_IN)
    gni=','.join(genomes_NOT_IN)
    if gni=='': 
        gni='None'
        not_in=False
    output.write('#ALL GENOMES\n#Genomes included :'+gi+' ; Genomes excluded :'+gni+'\n'+'#Parameters: PAM:'+reverse_complement(PAM)+' ; sgRNA size:'+str(non_PAM_motif_length)+'\n')
    output.write('sgRNA sequence')
    for genome_i in genomes_IN:
        output.write('\t'+genome_i) 
    output.write('\n')
    for hit in hit_list:
        output.write(hit.sequence)
        for gi in genomes_IN: 
            output.write('\t'+','.join(hit.genomes_Dict[gi]))              
        output.write('\n')  

def output_interface(hit_list,genomes_NOT_IN):   
    '''
    Reformat the results to print them in json format. 
    There will be parsed in javascript to display it in interface. 
    '''
    json='['
    for hit in hit_list: 
        json+='{"sequence":"'+hit.sequence+'","in":['
        for genome in hit.genomes_Dict:
            json+='{"org":"'+genome+'","coords":"'
            coords=';'.join(hit.genomes_Dict[genome])
            json+=coords+'"},'
        json=json.rstrip(',')    
        json+=']},'
    json=json.rstrip(',')    
    json+=']'    
    not_in_str=','.join(genomes_NOT_IN)  

    print(json)
    print(not_in_str)


def construction(indexs_path,fasta_path,bowtie_path,PAM,non_PAM_motif_length,genomes_IN,genomes_NOT_IN,dict_org_code):
    '''
    Launch the function needed that depends of users input, writes the results in file and creates a string 
    in json format for parsing the results in javascript.
    '''
    start_time=time.time()
    os.system('mkdir tmp')
    start = time.time()
    if len(genomes_IN)!=1:
        #hit_list=search_common_sgRNAs_by_construction(fasta_path,PAM,non_PAM_motif_length,genomes_IN,dict_org_code,bowtie_path,indexs_path)
        sorted_genomes=sort_genomes(genomes_IN,fasta_path,dict_org_code)
    else: 
        sorted_genomes=genomes_IN   

    dic_seq=construct_in(fasta_path,sorted_genomes[0],dict_org_code[sorted_genomes[0]],PAM,non_PAM_motif_length)
    eprint(str(len(dic_seq))+' hits in first included genome '+sorted_genomes[0])
    write_to_fasta(dic_seq)

    if len(genomes_NOT_IN)>=1: 
        sorted_genomes_notin=sort_genomes_desc(genomes_NOT_IN,fasta_path,dict_org_code)
        for genome in genomes_NOT_IN: 
            dic_seq=add_notin(dict_org_code[genome],indexs_path,dic_seq)
            if len(dic_seq)==0: 
                print("Program terminated&No hits remain after exclude genome "+genome)
                end_time=time.time()
                total_time=end_time-start_time
                eprint('TIME',total_time)
                sys.exit(1)
            eprint(str(len(dic_seq))+' hits remain after exclude genome '+genome) 
            write_to_fasta(dic_seq) 

    if len(sorted_genomes)>1: 
        for genome in sorted_genomes[1:]:
            dic_seq=add_in(bowtie_path,dict_org_code[genome],indexs_path,dic_seq,genome,len(PAM)+non_PAM_motif_length)
            if len(dic_seq)==0: 
                print("Program terminated&No hits remain after include genome "+genome)
                end_time=time.time()
                total_time=end_time-start_time
                eprint('TIME',total_time)
                sys.exit(1)
            eprint(str(len(dic_seq))+' hits remain after include genome '+genome)    
            write_to_fasta(dic_seq)

    hit_list=construct_hitlist(dic_seq)    

    hit_list=Scorage_triage(hit_list) 

    ##Put results in local file for access via the interface.
    write_to_file(genomes_IN,genomes_NOT_IN,hit_list[:1000],PAM,non_PAM_motif_length)

    ##Output formatting for printing to interface
    output_interface(hit_list[:100],genomes_NOT_IN)

    #os.system('rm -r tmp')


            
def main():
    start_time=time.time()
    indexs_path = './reference_genomes/index'
    fasta_path = './reference_genomes/fasta'
    bowtie_path='./bowtie-1.1.2/bowtie'
    dict_organism_code = construct_dict_organism_assemblyref()   ##Keys: organism, values: genomic reference (ncbi)
    dict_code_organism=intervert(dict_organism_code)      ##Keys/values interchanged relative to line above
    organisms_selected,organisms_excluded,PAM,non_PAM_motif_length=args_gestion(dict_organism_code)
    construction(indexs_path,fasta_path,bowtie_path,PAM,non_PAM_motif_length,organisms_selected,organisms_excluded,dict_organism_code)
    end_time=time.time()
    total_time=end_time-start_time
    eprint('TIME',total_time)
    eprint('CUMULATIVE_LENGTH',cumulative_length(organisms_selected,dict_organism_code))
    eprint('CUMULATIVE_LENGTH_NOT_IN',cumulative_length(organisms_excluded,dict_organism_code))
    
if __name__=='__main__':
    main()