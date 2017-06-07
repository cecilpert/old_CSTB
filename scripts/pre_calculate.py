from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
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
    print(reg_exp)
    indices = [m.start() for m in re.finditer('(?=' + reg_exp + ')', seq, re.I)]
    return indices

def find_sgRNA(organism_code,PAM,non_PAM_motif_length):
    fasta_file=fasta_path +'/' + organism_code +'_genomic.fna'
    genome_seqrecord=next(SeqIO.parse(fasta_file, 'fasta'))
    genome_seq=str(genome_seqrecord.seq)
    seq_list_forward=find_PAM(genome_seq,reverse_complement('NNNNNNNNNNNNNNNNNNNNNGG'))
    seq_list_reverse=find_PAM(genome_seq,'NNNNNNNNNNNNNNNNNNNNNGG')
    print(seq_list_forward,seq_list_reverse)
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
    return dic


def concatenate_seq(organism_code,dic_overlap): 
    fasta_file=fasta_path +'/' + organism_code +'_genomic.fna'
    genome_seqrecord=next(SeqIO.parse(fasta_file, 'fasta'))
    genome_seq=str(genome_seqrecord.seq)
    new_seq=''
    for start,end in dic_overlap.items(): 
        new_seq+=genome_seq[start:end]
    
    new_description=genome_seqrecord.description+',just sgRNA'
    new_seq_record=SeqRecord(Seq(new_seq),genome_seqrecord.id,genome_seqrecord.name,new_description)
    SeqIO.write(new_seq_record,'reference_genomes/pre_calculate/'+organism_code+'_sgRNA.fa','fasta')
      

fasta_path='./reference_genomes/fasta'
dic_overlap=find_sgRNA('GCF_000007325.1_ASM732v1','NGG',20)
concatenate_seq('GCF_000007325.1_ASM732v1',dic_overlap)