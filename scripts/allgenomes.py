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
        self.NOTIN_Dict={}
        self.score=()   ##Holds score as tuple of occurence score (on genomes_IN) and NOTIN score 


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


def construct_seq_dict(fasta_path,organism_code,PAM,non_PAM_motif_length):
    '''
    Construct seq dict for an organism. 
    1. Take the nt sequence for the genome as SeqRecord object
    2. Search PAM motif (for strand +) and reverse complement PAM motif (for strand -) in the sequence with find_sgRNA_seq function. This function give a list with indexes where the motif was found.
    4. Browse index lists and take the sequence of the complete sgRNA from index. 
    To determinate the sequence of the complete sgRNA : 
    - check if the motif isn't to close from start or end of the genome sequence 
    - determinate coordinates, see code for details 
    - take the subsequence in genome SeqRecord with calculate coordinates. It will be reverse complemented for strand +. 
    - create dictionnary like key=sequence and value=list of coordinates like strand(start,end)
    5. Return the dictionnary
    '''

    fasta_file=fasta_path +'/' + organism_code +'_genomic.fna'
    genome_seqrecord=next(SeqIO.parse(fasta_file, 'fasta'))

    genome_length=len(genome_seqrecord)

    genome_seq=str(genome_seqrecord.seq)
    
    seq_list_forward=find_sgRNA_seq(genome_seq,PAM)
    seq_list_reverse=find_sgRNA_seq(genome_seq,reverse_complement(PAM))

    seq_dict={}

    count=0
    for indice in seq_list_forward: 
        end=indice+len(PAM)+non_PAM_motif_length
        if not end > genome_length:
            seq=str(genome_seqrecord.seq[indice:end].reverse_complement())
            if seq not in seq_dict:
                count+=1
                seq_dict[seq]=['+('+str(indice+1)+','+str(end)+')']
            else:   
                seq_dict[seq].append('+('+str(indice+1)+','+str(end)+')')
            
    for indice in seq_list_reverse: 
        end=indice+len(PAM)
        start=indice-non_PAM_motif_length
        if not start<=0:
            seq=genome_seq[start:end]
            if seq not in seq_dict:
                seq_dict[seq]=['-('+str(start+1)+','+str(end)+')'] 
            else:    
                seq_dict[seq].append('-('+str(start+1)+','+str(end)+')')
          
    return(seq_dict)    


def construct_seq_dict_with_organism(fasta_path,organism,organism_code,PAM,non_PAM_motif_length):
    '''
    Same function as construct_seq_dict except the dictionnary will not be like value=list of coordinates. It add the information about organism with dictionnary values like : dictionnary with key=organism and value=list of coordinates in this organism
    '''
    start=time.time()
    fasta_file=fasta_path +'/' + organism_code +'_genomic.fna'
    genome_seqrecord=next(SeqIO.parse(fasta_file, 'fasta'))

    genome_length=len(genome_seqrecord)
    eprint('length',genome_length)

    genome_seq=str(genome_seqrecord.seq)
    
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

   
def search_common_sgRNAs_by_construction(fasta_path,PAM,non_PAM_motif_length,genomes_IN,dict_org_code,bowtie_path,indexs_path):
    '''Create dictionnary with common sgRNA between all organisms and call the function to create hitlist for this dictionnary. The final dictionnary is like key=sequence and value=other dictionnary with key=organism and value=list of position in the organism
    Steps : 
    1. Take the smallest genome for the first iteration because we can't have more common sgRNA than total sgRNA of the smallest genome 
    2. Create a dictionnary of smallest genome like key=sequence and value=dictionnary with key=organism and value=list of position
    3. Write all sequences in dictionnary in a temporary fasta file
    4. Use bowtie with this file as query and the second genome as reference 
    5. Use the function treat_bowtie_results to generate a dictionnary with key=sequence and value=position in genome.
    6. Complete the dictionnary the first dictionnary that included informations about organism by adding the position in the second genomes. 
    6. Write the common sequences in temporary fasta file and launch bowtie on next genome. Do this steps until all included genomes are treated. 
    7. Use the function construct_hit_list that will construct an hit list with final dictionnary.  
    '''

    out=open('tmp/sgRNA_in.fa','w')
    new_genomes_IN=[]

    for element in genomes_IN:
        fasta_genome=next(SeqIO.parse(fasta_path +'/' + dict_org_code[element] +'_genomic.fna', 'fasta'))
        new_genomes_IN.append((len(fasta_genome.seq),element)) ##Tuple with length of genome and genome reference code

    genomes_IN=[i[1] for i in sorted(new_genomes_IN,key=lambda genome:genome[0])]   ##Sort by ascending size
    smallest=genomes_IN[0]
    #eprint("Sorted genomes by size, extracted smallest Dict,",smallest)
    ##IMPROVEMENT: SORT GENOMES_IN BY INCREASING SGRNA DICT SIZE

    start_global=time.time()
    start=time.time()
    seq_dict_to_compare= construct_seq_dict_with_organism(fasta_path,smallest,dict_org_code[smallest],PAM,non_PAM_motif_length)
    end=time.time()
    print('first_hits',len(seq_dict_to_compare))
    construction_time=end-start

    for seq in seq_dict_to_compare: 
        out.write('>'+seq+'\n'+seq+'\n')

    out.close()    

    for organism in genomes_IN[1:]:
        new_dic={}
        organism_code=dict_org_code[organism]
        bowtie_command=bowtie_path+' --quiet -a -c -v 0 '+indexs_path+'1/'+organism_code+' -f tmp/sgRNA_in.fa tmp/results_bowtie.txt'
        os.system(bowtie_command)
        dic_results=treat_bowtie_results_IN(PAM,non_PAM_motif_length)
        out=open('tmp/sgRNA_in.fa','w')
        for seq in dic_results: 
            out.write('>'+seq+'\n'+seq+'\n')
            new_dic[seq]=seq_dict_to_compare[seq]
            new_dic[seq][organism]=dic_results[seq]    
        out.close()
        seq_dict_to_compare=new_dic       
        if not seq_dict_to_compare:  # Empty list evaluates to False
            print("Program terminated& Adding organism "+organism+" gave an absence of common occurences. Maybe try lowering sgRNA size. \nCAVEAT: the search is for exact matches.")
            end_global=time.time()  
            time_global=end_global-start_global
            print('TIME',time_global)
            sys.exit(1)
    end_global=time.time()  
    time_global=end_global-start_global   
    eprint(time_global)   
    print('CONSTRUCTION',construction_time)
    print('COMPARAISON',end_global-start_global-construction_time)
    hits_list=construct_hitlist(seq_dict_to_compare)

    return hits_list

def treat_bowtie_results_IN(PAM,non_PAM_motif_length): 
    '''Treat bowtie results for included genomes. 

    Take submit sequences, strand and pos of alignment. 
    Store informations in dictionnary with key=sequence and value=list of coordinates in genome (format of coordinates is : strand(start,end))
    Return this dictionnary 
    '''
    f=open('tmp/results_bowtie.txt','r')
    res_dict={}    
    len_sgRNA=len(PAM)+int(non_PAM_motif_length)
    forbidden_positions=[i for i in range(int(non_PAM_motif_length),len_sgRNA)]
    dic_strand={'+':'-','-':'+'}
    dic_results={}
    for line in f: 
        line_split=line.rstrip().split('\t')
        hit_seq=line_split[0]   
        strand=line_split[1]
        pos=int(line_split[3])
        if strand=='+': 
            start=pos+1
            end=pos+len_sgRNA
        elif strand=='-': 
            end=pos
            start=pos-len_sgRNA
        if hit_seq not in dic_results:     
            dic_results[hit_seq]=[dic_strand[strand]+'('+str(start)+','+str(end)+')']  
        else: 
            dic_results[hit_seq].append(dic_strand[strand]+'('+str(start)+','+str(end)+')')
    f.close()
    return dic_results 

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


def not_in_search(indexs_path,bowtie_path,genomes_NOT_IN,dict_org_code,hit_list,PAM,non_PAM_motif_length):
    '''Search against excluded genomes.
    - For each genomes, launch bowtie2 with all common sgRNA as query and the genome as reference. The length of the seed (option -L) is changed because the default size (22) is too close of the size of sgRNA to allow enough mismatches in research. The option -a write all alignments in bowtie results.
    - Bowtie results are treated with treat_bowtie_results_not_in function. This function will complete each sgRNA object Hit with information about alignments with excluded genomes.
    '''
    for organism in genomes_NOT_IN:
        organism_code=dict_org_code[organism]
        bowtie_command='bowtie2 -x'+indexs_path+'2/'+organism_code+' -f tmp/sgRNA_in.fa -S tmp/results_bowtie2.sam -L 13 -a'
        os.system(bowtie_command)
        treat_bowtie_results_not_in(PAM,organism,hit_list)

def treat_bowtie_results_not_in(PAM,organism,hit_list): 
    '''
    Treat bowtie2 results for excluded genomes. 
    Split results and take submits sequences, sequences that align, cigar code and position of alignment. 
    Treat this informations and store the results in a dictionnary (key=sequence and value=information about alignment). Information about alignment will be like : mismatch_information:strand and position:supplemetary information. The first field can be 'Exact match', 'No indel' (it means no indel but mismatch) or 'Indel' with information about indel position. The third field contains informations about mismatches positions. 
    Example of indel position (cigar code) : 15M1I7M. This means that there is 15 aligned position (match or mismatch), 1 insertion in submit sequence and 7 aligned positions.
    Example of mismatch position : 15A8. This means the 15 first base are the same, then there is a A in submit sequence and an other base in reference sequence, and the 8 last bases are the same again. 
    See more about cigar code at : https://samtools.github.io/hts-specs/SAMv1.pdf (Section 1.4) 
    '''
    f=open('tmp/results_bowtie2.sam','r')
    res_dict={}
    for l in f: 
        if l[0]!='@':
            l_split=l.rstrip().split('\t')
            hit_seq=l_split[0]
            seq_align=l_split[9]
            cigar=l_split[5]
            pos=l_split[3]
            if cigar=='*': 
                res_dict[hit_seq]='None'
            else: 
                if hit_seq not in res_dict: 
                    res_dict[hit_seq]=[]
                mm=l_split[-2]
                mm=mm.split(':')[-1]
                if seq_align==hit_seq: 
                    strand='-'
                elif seq_align!=hit_seq:
                    strand='+' 
                    cigar=reverse_cigar(cigar)   
                if mm=='23' and cigar=='23M':
                    res_dict[hit_seq].append('Exact match:'+strand+pos)
                else: 
                    if cigar=='23M': 
                        dic_str='No indel:'+strand+pos
                    else: 
                        dic_str='Indel '+cigar+':'+strand+pos    
                    count_mm=0
                    if strand=='-': 
                        mm_pam=treat_mismatch_reverse(PAM,mm) 
                    elif strand=='+':  
                        mm_pam=treat_mismatch_forward(PAM,mm)  
                        mm=reverse_mismatch(mm)                         
                    if mm_pam: 
                        res_dict[hit_seq].append(dic_str+':Mismatch in PAM')
                    else: 
                        res_dict[hit_seq].append(dic_str+':'+mm)
                                                                                                                        
    f.close()             
    complete_hit_list(hit_list,res_dict,organism)      

def reverse_mismatch(mm): 
    '''
    Take a mismatch string as input and return its reverse. The string is read from the end and the base are complemented. 
    Example : for 10A2A9C2, the reverse is 2G9T2T10
    '''
    new_mm=''
    rc_dic={'A':'T','C':'G','G':'C','T':'A','N':'N'}
    list_bases=[]
    for i in mm: 
        if i.isalpha(): 
            list_bases.append(i)   
    list_bases.reverse()          
    for base in list_bases: 
        mm_split=mm.split(base)
        new_mm+=mm_split[-1]+rc_dic[base]
        mm=mm_split[:-1]
        if len(mm)>1: 
            mm=base.join(mm)
        else:
            mm=mm[0]   
        if mm.isnumeric():
            new_mm+=mm
    return(new_mm)     

def reverse_cigar(cigar):
    '''
    Reverse a cigar code. 
    Example : for 13M1I9M, the reverse is 9M1I13M
    '''
    new_cigar=''
    list_letters=[]
    for i in cigar: 
        if i.isalpha(): 
            list_letters.append(i)
    first=list_letters.pop()
    list_letters.reverse()
    for letter in list_letters: 
        cigar_split=cigar.split(letter)
        new_cigar+=cigar_split[-1]
        cigar=cigar_split[:-1]         
        if len(cigar)>1: 
            cigar=letter.join(cigar)+letter
        else:
            cigar=cigar[0]+letter   
    new_cigar+=cigar
    return new_cigar        

           
def treat_mismatch_forward(PAM,mm): 
    '''
    Check if a mismatch is in the PAM motif. For forward verification, check if the mismatch is in the first bases (depending on motif length).  
    '''
    acgt=['A','C','G','T']
    for base in acgt: 
        mm_pam=mm.split(base)[0]
        if mm_pam.isnumeric():
            if int(mm_pam)<len(PAM):
                return True
    return False        

def treat_mismatch_reverse(PAM,mm):
    '''
    Check if a mismatch is in the PAM motif. For reverse verification, check if the mismatch is in the last bases (depending on motif length).  
    '''
    acgt=['A','C','G','T']
    for base in acgt: 
        mm_pam=mm.split(base)[-1]
        if mm_pam.isnumeric():
            if int(mm_pam)<len(PAM):
                return True
    return False  
                   
                          
def complete_hit_list(hit_list,res_dict,organism): 
    '''Complete Hit objects previously created with construct_hit_list function.
    Add informations about positions of the sequences in excluded genomes '''
    for hit in hit_list: 
        if hit.sequence in res_dict: 
            hit.NOTIN_Dict[organism]=res_dict[hit.sequence]       

def Scorage_triage(hitlist):
    '''
    Scoring of the hits found, where positive scores mean stronger sgRNA constructs.
    Complete Hit objects with score and sort the hits with this scores. 
    '''
    parametre_triage=100
    sorted_hitlist=[]
    NOTIN_absence=1000   ##Scorage definitions. Absolute priority on NOTIN absence followed by sorting through exactmatch/mismatches.
    NOTIN_sub=0.1
    IN_occ=1
    for hit in hitlist:
        IN_score=0
        NOTIN_score=0
        for j in hit.genomes_Dict:
            IN_score+=(len(hit.genomes_Dict[j])*IN_occ)
        for j in hit.NOTIN_Dict:    ##Will not be evaluated if no NOTIN search as NOTIN_Dict will be empty.
            if hit.NOTIN_Dict[j]=="None":
                NOTIN_score+=NOTIN_absence
            else:
                count=0
                for k in hit.NOTIN_Dict[j]:
                    count-=1    ##The match of an sgRNA on a NOTIN org is down-scored.
                    if k.split(':')[0]!='Exact match':   ##Implicitly, an exact match is given a score of 0.
                        NOTIN_score+=number_mismatch(k)*NOTIN_sub
                NOTIN_score+=count
        hit.score=(IN_score,NOTIN_score)
    sorted_hitlist=sorted(hitlist,key=lambda hit:(hit.score[1]+hit.score[0])/2,reverse=True)
    return(sorted_hitlist)

def number_mismatch(mm): 
    if mm=='-Mismatch in PAM' or mm=='+Mismatch in PAM': 
        return 0.1
    else: 
        count_mm=0
        mm=mm.lstrip('-').lstrip('+')  
        for i in mm: 
            if i.isalpha(): 
                count_mm+=1        
        return(count_mm)   
             
def write_to_file(genomes_IN,genomes_NOT_IN,hit_list,PAM,non_PAM_motif_length): 
    '''Write results in a file. 
    The file is a tabulated file, with first column=sgRNA sequence, then one column for each included organisms with list of positions of the sequence in each, and finally one column for each excluded organism with informations about alignment of the sequence with this genomes organisms.
    '''
    new_tag=timestamp()
    output=open('./scripts/static/results'+new_tag+'.txt','w')
    not_in=True
    gi=','.join(genomes_IN)
    gni=','.join(genomes_NOT_IN)
    if gni=='': 
        gni='None'
        not_in=False
    print(gi)
    print(gni)
    output.write('#ALL GENOMES\n#Genomes included :'+gi+' ; Genomes excluded :'+gni+'\n'+'#Parameters: PAM:'+reverse_complement(PAM)+' ; sgRNA size:'+str(non_PAM_motif_length)+'\n')
    output.write('sgRNA sequence')
    for genome_i in genomes_IN:
        output.write('\t'+genome_i)
    if not_in: 
        for genome_ni in genomes_NOT_IN: 
            output.write('\t'+genome_ni)   
    output.write('\n')
    for hit in hit_list:
        output.write(hit.sequence)
        for gi in genomes_IN: 
            output.write('\t'+','.join(hit.genomes_Dict[gi]))
        if not_in: 
            for gni in genomes_NOT_IN: 
                if hit.NOTIN_Dict[gni]=='None':
                    output.write('\t Absent')
                else: 
                    output.write('\t')
                    exact_match_counter=0
                    mismatch_arr=[]
                    exact_match_arr=[]
                    for el in hit.NOTIN_Dict[gni]:
                        if el.split(':')[0]=='Exact match':
                            exact_match_arr.append(el.split(':')[1])
                        else:
                            mismatch_arr.append(el)
                    if len(exact_match_arr)>0 and len(mismatch_arr)>0: 
                        to_write=str(len(exact_match_arr))+' exact match(es) : '+','.join(exact_match_arr)+ ' and '+str(len(mismatch_arr))+' non-exact match(es): '+','.join(mismatch_arr)    
                    elif len(exact_match_arr)>0: 
                       to_write=str(len(exact_match_arr))+' exact match(es) : '+','.join(exact_match_arr)
                    elif len(mismatch_arr)>0: 
                        to_write=str(len(mismatch_arr))+' non-exact match(es): '+','.join(mismatch_arr)   
                    output.write(to_write)               
        output.write('\n')  


def output_interface(hit_list,genomes_NOT_IN,new_tag):   
    '''
    Reformat the results to print them in json format. 
    There will be parsed in javascript to display it in interface. 
    '''
    output_arr=[]
    output_arr_NOTIN=[]
    seq_arr=[]
    org_arr=[]
    org_arr_NOTIN=[]
    counter=0
    for i in hit_list: 
        genome_arr=[]
        genome_arr_NOTIN=[]
        seq_arr.append(i.sequence) 

        ##Hybridisation output
        for j in i.genomes_Dict:
            if counter==0:
                org_arr.append(j)
            occ_arr=[]
            for k in range(len(i.genomes_Dict[j])):
                occ_arr.append(str(i.genomes_Dict[j][k]))
            genome_arr.append('{"'+j+'":"'+'%'.join(occ_arr)+'"}')

        ##NOT IN output
        for NOTIN in i.NOTIN_Dict:
            if i.NOTIN_Dict[NOTIN]=="None":
                occ_NOTIN='Absent'
            else:
                occ_NOTIN=''
                occ_NOTIN_arr=i.NOTIN_Dict[NOTIN]
                mismatch_arr=[]
                exact_match_arr=[]
                for el in occ_NOTIN_arr:
                    if el.split(':')[0]=='Exact match':
                        exact_match_arr.append(el)
                    else:
                        mismatch_arr.append(el)
                if mismatch_arr:
                    occ_NOTIN+='Present:'+str(len(exact_match_arr))+' exact match(es):'+','.join(exact_match_arr) +' and '+str(len(mismatch_arr))+' non-exact matches: '+','.join(mismatch_arr)
                else:
                    occ_NOTIN+='Present:'+str(len(exact_match_arr))+' exact match(es):'+','.join(exact_match_arr)
            genome_arr_NOTIN.append('{"'+NOTIN+'":"'+occ_NOTIN+'"}') 

        output_arr.append('{"'+i.sequence+'":['+','.join(genome_arr)+']}')
        output_arr_NOTIN.append('{"'+i.sequence+'":['+','.join(genome_arr_NOTIN)+']}')
        counter+=1
    for nig in genomes_NOT_IN:
        org_arr_NOTIN.append(nig)
    result_string='{"sequences":['+','.join(output_arr)+']}'
    result_string_NOTIN='{"sequences":['+','.join(output_arr_NOTIN)+']}'
    seq_string=','.join(seq_arr)
    org_string=','.join(org_arr)
    if org_arr_NOTIN:
        org_NOTIN_string=','.join(org_arr_NOTIN)
    else:
        org_NOTIN_string='None'
    print(new_tag)
    print(org_string)
    print(org_NOTIN_string)
    print(result_string)
    print(result_string_NOTIN)
    print(seq_string) 
              
def construction(indexs_path,fasta_path,bowtie_path,PAM,non_PAM_motif_length,genomes_IN,genomes_NOT_IN,dict_org_code):
    '''
    Launch the function needed that depends of users input, writes the results in file and creates a string 
    in json format for parsing the results in javascript.
    '''
    os.system('mkdir tmp')
    start = time.time()
    if len(genomes_IN)!=1:
        hit_list=search_common_sgRNAs_by_construction(fasta_path,PAM,non_PAM_motif_length,genomes_IN,dict_org_code,bowtie_path,indexs_path)
    else:
        start_const=time.time()
        organism=genomes_IN[0]
        seq_dict=construct_seq_dict_with_organism(fasta_path,organism,dict_org_code[organism],PAM,non_PAM_motif_length)
        hit_list = construct_hitlist(seq_dict)
        end_const=time.time()
    eprint('Found',len(hit_list),'hits on',len(genomes_IN),'genomes using an sgRNA of size',
          len(PAM)+non_PAM_motif_length,'in',time.time()-start,'seconds.')
    not_in_status="NOT IN status: "
    if len(genomes_NOT_IN)==0:
        not_in_status+="No"
        eprint(not_in_status)
    else:
        not_in_status+="Yes"
        eprint(not_in_status)
        not_in_search(indexs_path,bowtie_path,genomes_NOT_IN,dict_org_code,hit_list,PAM,non_PAM_motif_length)
        
    #Score and sort the results       
    hit_list=Scorage_triage(hit_list) 

    ##Put results in local file for access via the interface.
    new_tag=write_to_file(genomes_IN,genomes_NOT_IN,hit_list[:1000],PAM,non_PAM_motif_length)

    ##Output formatting for printing to interface
    output_interface(hit_list[:100],genomes_NOT_IN,new_tag)

    os.system('rm -r tmp')


            
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