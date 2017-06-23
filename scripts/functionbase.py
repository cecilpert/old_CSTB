from __future__ import print_function
import time,os
from Bio import SeqIO

def construct_dict_organism_assemblyref():
    try:
        genomes = open('./scripts/reference_genomes_ftp.txt','r')
    except:
        print("File opening error")
        exit(1)

    dic = {}
    for genome in genomes:
        genome = genome.rstrip()
        columns = genome.split("\t")
        if columns[0] not in dic:
            dic[columns[0]] = columns[1]

    genomes.close()
    org_assembly_Dict = {}
    for i in dic:
        org_assembly_Dict[i] = dic[i].split('/')[-1]
    return org_assembly_Dict

def construct_list_organism_ftp():
    try:
        genomes = open('./scripts/reference_genomes_ftp.txt','r')
    except:
        print("File opening error")
        exit(1)
    org_assembly_list=[]
    for genome in genomes:
        genome = genome.rstrip()
        columns = genome.split("\t")
        assembly_reference=columns[1].split('/')[-1]
        org_assembly_list.append([columns[0],assembly_reference])
    org_assembly_list.sort(key=lambda x: x[0])
    return org_assembly_list

def intervert(Dict):
    newDict={}
    for i in Dict:
        newDict[Dict[i]]=i
    return newDict

def reverse_complement(sequence):
    '''
    Function for turning a 5'-3' nucleotidic sequence into its 5'-3' reverse complement.
    '''
    rev_comp = []
    for idx in range(len(sequence) - 1, -1, -1):
        if sequence[idx] == 'A':
            rev_comp = rev_comp + ['T']
        elif sequence[idx] == 'C':
            rev_comp = rev_comp + ['G']
        elif sequence[idx] == 'G':
            rev_comp = rev_comp + ['C']
        elif sequence[idx] == 'T':
            rev_comp = rev_comp + ['A']
        else:
            rev_comp = rev_comp + ['N']
    return "".join(rev_comp)

def build_expression(seq):
    result = ''
    iupac_code={'R':'[AG]', 'Y':'[CT]', 'S':'[GC]', 'W':'[AT]', 'K':'[GT]', 'M':'[AC]', 'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]'}
    for c in seq:
        if c in iupac_code:
            result = result + iupac_code[c]
        else:
            result = result + c
    return result

def timestamp():
    '''
    Generates a number corresponding to the current time (a float representation of the time since year 1970), used for appending to the results file so that it will not
    have been cached by the browser. Note the current time is always increasing, EXCEPT IF YOU RESET YOUR SYSTEM CLOCK.
    '''
    tag=open('./scripts/static/resultsTag.txt','r')
    content=tag.readline().rstrip()
    tag.close()
    os.remove('./scripts/static/results'+content+'.txt')
    new_tag=str(time.time())
    tag=open('./scripts/static/resultsTag.txt','w')
    tag.write(new_tag)
    tag.close()
    return(new_tag)

def reverse_complement_mismatch(i,len_sgrna):
    '''Reverse complement mismatch to give position of mismatch on sequence NNNNN+PAM (sequence printed in results). Mismatch is given in 1-based indexing.
    "i" is the original mismatch, in '7:A>G' format.
    '''
    n_i=i.split(':')
    n_i[0]=str(len_sgrna-int(n_i[0])+1)
    n_i[1]='>'.join([reverse_complement(j) for j in n_i[1].split('>')])
    new_mismatch=':'.join(n_i)
    return new_mismatch

def corrected_mismatch(mm,len_sgrna,strand):
    mm_split=mm.split(':')
    mm_split[0]=str(len_sgrna-int(mm_split[0]))
    if strand=='+':
        mm_split[1]='>'.join([reverse_complement(j) for j in mm_split[1].split('>')])
    else:
        mm_split[1]='>'.join([j for j in mm_split[1].split('>')])
    new_mismatch=':'.join(mm_split)
    return new_mismatch

def cumulative_length(genomes_in,dict_org_code):
    '''Function used for profiling the code, return the cumulative length of all genomes given'''
    cumulative_length=0
    for genome in genomes_in:
        count=0
        code=dict_org_code[genome]
        for seq_record in SeqIO.parse('reference_genomes/fasta/'+code+'_genomic.fna','fasta'):
            if count==0:
                cumulative_length=cumulative_length+len(seq_record)
            count+=1
    return cumulative_length