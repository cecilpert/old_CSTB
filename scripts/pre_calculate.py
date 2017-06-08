from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re, pickle, os
from functionbase import reverse_complement

def build_expression(seq):
    result = ''
    iupac_code={'R':'[AG]', 'Y':'[CT]', 'S':'[GC]', 'W':'[AT]', 'K':'[GT]', 'M':'[AC]', 'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]'}
    for c in seq:
        if c in iupac_code:
            result = result + iupac_code[c]
        else:
            result = result + c
    return result

def find_PAM(seq,motif): 
    list_seq=[]
    reg_exp = build_expression(motif) 
    indices = [m.start() for m in re.finditer('(?=' + reg_exp + ')', seq, re.I)]
    return indices

def find_sgRNA(organism_code,PAM,non_PAM_motif_length):
    fasta_file='reference_genomes/fasta/' + organism_code +'_genomic.fna'
    genome_seqrecord=next(SeqIO.parse(fasta_file, 'fasta'))
    genome_seq=str(genome_seqrecord.seq)
    sgRNA='' 
    for i in range(non_PAM_motif_length): 
        sgRNA+='N'
    sgRNA+=PAM    
    seq_list_forward=find_PAM(genome_seq,reverse_complement(sgRNA))
    seq_list_reverse=find_PAM(genome_seq,sgRNA)

    return seq_list_forward,seq_list_reverse

def concatenate_seq(seq_list_forward,seq_list_reverse,organism_code):    
    seq_list_all=seq_list_forward+seq_list_reverse
    seq_list_all=set(seq_list_all)
    seq_list_all=list(seq_list_all)
    seq_list_all.sort()

    previous_end=0
    previous_start=0
    count_overlap=0
    previous_start_to_keep=seq_list_reverse[0]
    dic={}
    for i in seq_list_reverse: 
        start=i 
        end=i+23
        if start<previous_end: 
            count_overlap+=1
            overlap=True
            if count_overlap==1: 
                start_to_keep=previous_start 
        else: 
            count_overlap=0 
            overlap=False
        if not overlap: 
            start_to_keep=start
        end_to_keep=end     
        dic[start_to_keep]=end_to_keep
        previous_end=end  
        previous_start=start  
    

    fasta_file='reference_genomes/fasta/' + organism_code +'_genomic.fna'
    genome_seqrecord=next(SeqIO.parse(fasta_file, 'fasta'))
    genome_seq=str(genome_seqrecord.seq)
    new_seq=''
    for start,end in dic_overlap.items(): 
        new_seq+=genome_seq[start:end]
    
    new_description=genome_seqrecord.description+',just sgRNA'
    new_seq_record=SeqRecord(Seq(new_seq),genome_seqrecord.id,genome_seqrecord.name,new_description)
    SeqIO.write(new_seq_record,'reference_genomes/pre_calculate/'+organism_code+'_sgRNA.fa','fasta')

def store_positions(seq_list_forward,seq_list_reverse,organism_code,PAM,non_PAM_motif_length):
    '''
    Same function as construct_seq_dict except the dictionnary will not be like value=list of coordinates. It add the information about organism with dictionnary values like : dictionnary with key=organism and value=list of coordinates in this organism
    '''
    fasta_file='reference_genomes/fasta/' + organism_code +'_genomic.fna'
    genome_seqrecord=next(SeqIO.parse(fasta_file, 'fasta'))
    genome_seq=str(genome_seqrecord.seq)
    seq_dict={}
    for indice in seq_list_forward: 
        end=indice+len(PAM)+non_PAM_motif_length
        seq=str(genome_seqrecord.seq[indice:end].reverse_complement())
        if seq not in seq_dict:
            seq_dict[seq]=[]
        seq_dict[seq].append('+('+str(indice+1)+','+str(end)+')')
            
    for indice in seq_list_reverse: 
        end=indice+len(PAM)
        start=indice-non_PAM_motif_length
        seq=genome_seq[start:end]
        if seq not in seq_dict:
            seq_dict[seq]=[]   
        seq_dict[seq].append('-('+str(start+1)+','+str(end)+')')
    
    pickle.dump(seq_dict,open("reference_genomes/pre_calculate/"+organism_code+"_dicpos.pic","wb"))

    to_write=''
    for seq in seq_dict: 
        to_write+=seq+'\t'+';'.join(seq_dict[seq])+'\n'

    out=open("reference_genomes/pre_calculate/"+organism_code+"_pos.txt",'w')
    out.write(to_write)
    out.close()

    
      

def launch_all_genomes(file_all_genomes,PAM,non_PAM_motif_length): 
    f=open(file_all_genomes,'r')     
    for l in f: 
        organism_code=l.rstrip().split('\t')[1].split('/')[-1]
        print(organism_code)
        list_forward,list_reverse=find_sgRNA(organism_code,PAM,non_PAM_motif_length)
        concatenate_seq(list_forward,list_reverse,organism_code)


#launch_all_genomes('scripts/reference_genomes_ftp.txt','NGG',20)
list_test=['GCF_000005845.2_ASM584v2','GCF_000008865.1_ASM886v1','GCF_000026325.1_ASM2632v1','GCF_000026345.1_ASM2634v1','GCF_000183345.1_ASM18334v1','GCF_000299455.1_ASM29945v1','GCF_000012005.1_ASM1200v1','GCF_000006925.2_ASM692v2','GCF_000240185.1_ASM24018v2','GCF_000006945.1_ASM694v1','GCF_000195995.1_ASM19599v1','GCF_000215745.1_ASM21574v1','GCF_000025565.1_ASM2556v1']
for organism_code in list_test: 
    #print(organism_code)
    #list_forward,list_reverse=find_sgRNA(organism_code,'NGG',20)
    #store_positions(list_forward,list_reverse,organism_code,'NGG',20)
    index_cmd='bowtie2-build reference_genomes/pre_calculate/'+organism_code+'_sgRNA.fa reference_genomes/pre_calculate/'+organism_code+'_sgRNA'
    os.system(index_cmd)