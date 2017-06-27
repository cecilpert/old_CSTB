#from __future__ import print_function
from Bio.Seq import Seq
from Bio import SeqIO
from functionbase import *
from threading import Thread
from multiprocessing import Process
from queue import Queue
import time,argparse,os,sys,re,random
import subprocess
import pickle
import json

import uuid

#Global variable
TASK_KEY = str(uuid.uuid1())
WORKDIR = None
REF_GEN_DIR = None

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)     ##Found at http://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python


class Hit():
    def __init__(self, sequence, genomes_Dict):
        self.sequence=sequence
        self.genomes_Dict=genomes_Dict
        self.score=0  


def args_gestion():
    '''Take and treat arguments that user gives in command line'''
    ##Argparsing
    parser=argparse.ArgumentParser(description="Allgenomes program.")
    parser.add_argument("-gi",metavar="<str>",help="The organisms to search inclusion in.",required=True)
    parser.add_argument("-gni",metavar="<str>",help="The organisms to search exclusion from",required=True)
    parser.add_argument("-sl",metavar="<int>",help="The length of the sgRNA, excluding PAM")
    parser.add_argument("-pam",metavar="<str>",help="The PAM motif",required=True)
    parser.add_argument("-rfg",metavar="<str>",help="The path to the reference genome folder")
    parser.add_argument("-cah",metavar="<str>",help="The path to cache folder for the webservice")
    args=parser.parse_args()
    ##
    return args


def setupApplication(parameters, dict_organism_code):
    '''Processes the arguments to be usable for research'''
    organisms_selected=parameters.gi.split('+')
    organisms_excluded=parameters.gni.split('+')
    organisms_selected=[i for i in organisms_selected if i in dict_organism_code]
    organisms_excluded=[i for i in organisms_excluded if i in dict_organism_code]
    non_PAM_motif_length=int(parameters.sl)

    return organisms_selected, organisms_excluded, parameters.pam, non_PAM_motif_length

def setupWorkSpace(parameters):
    '''Create the work folder where all temporary or results files will be stored'''
    workFolder = parameters.cah + '/' + TASK_KEY
    os.mkdir(workFolder)
    return workFolder


def readOrganismCode(path):
    '''Load the dictionnary that contains genomes name associated with their references'''
    with open(path) as json_data:
        d = json.load(json_data)
    return d

def find_sgRNA_seq(seq,pam):
    '''
    Uses Regular expression matching of the PAM motif to the reference genome to get the start
    positions (0-based) of each match
    '''
    list_seq=[]
    reg_exp = build_expression(pam)
    indices = [m.start() for m in re.finditer('(?=' + reg_exp + ')', seq, re.I)]
    return indices


def construct_in(fasta_path, organism, organism_code, PAM, non_PAM_motif_length):
    '''Construct the sequences for first organism, with python regular expression research'''
    fasta_file=fasta_path + '/' + organism_code + '/' + organism_code +'_genomic.fna'

    sgRNA=''
    for i in range(non_PAM_motif_length):
        sgRNA+='N'
    sgRNA+=PAM

    seq_dict={}

    for genome_seqrecord in SeqIO.parse(fasta_file,'fasta'): 
        genome_seq=genome_seqrecord.seq
        ref=genome_seqrecord.id
        seq_list_forward=find_sgRNA_seq(str(genome_seq),reverse_complement(sgRNA))
        seq_list_reverse=find_sgRNA_seq(str(genome_seq),sgRNA)

        for indice in seq_list_forward:
            end=indice+len(PAM)+non_PAM_motif_length
            seq=genome_seq[indice:end].reverse_complement()
            seq=str(seq)
            if seq not in seq_dict:
                seq_dict[seq]={organism:[]}
            seq_dict[seq][organism].append(ref+':+('+str(indice+1)+','+str(end)+')')

        for indice in seq_list_reverse:
            end=indice+len(PAM)+non_PAM_motif_length
            seq=genome_seq[indice:end]
            seq=str(seq)
            if seq not in seq_dict:
                seq_dict[seq]={organism:[]}
            seq_dict[seq][organism].append(ref+':-('+str(indice+1)+','+str(end)+')')

    return seq_dict


def sort_genomes(list_genomes,fasta_path,dict_org_code):
    '''Sort genomes by ascending size'''
    tmp_list=[]
    for genome in list_genomes:
        len_genome=0
        for seq_record in SeqIO.parse(fasta_path +'/' + dict_org_code[genome][0] +'/' + dict_org_code[genome][0] +'_genomic.fna', 'fasta'): 
            len_genome+=len(seq_record)
        tmp_list.append((len_genome,genome))
    genomes_sorted=[i[1] for i in sorted(tmp_list,key=lambda genome:genome[0])]   ##Sort by ascending size
    return(genomes_sorted)

def sort_genomes_desc(list_genomes,fasta_path,dict_org_code):
    '''Sort genomes by descending size'''
    tmp_list=[]
    for genome in list_genomes:
        len_genome=0
        for seq_record in SeqIO.parse(fasta_path +'/' + dict_org_code[genome][0] +'/' + dict_org_code[genome][0]+'_genomic.fna', 'fasta'): 
            len_genome+=len(seq_record)
        tmp_list.append((len_genome,genome))
    genomes_sorted=[i[1] for i in sorted(tmp_list,key=lambda genome:genome[0],reverse=True)]   ##Sort by ascending size
    return(genomes_sorted)

def construct_hitlist(dict_seq):
    '''
    Will construct an object Hit for each sequence in dictionnary, and store all Hits in a list.
    These function only fill the attributes sequence and genomes_Dict of the object
    '''
    eprint('-- Construct final list --')
    hits_list=[]
    count=0
    for seq in dict_seq:
        count+=1
        new_hit=Hit(seq,dict_seq[seq])
        hits_list.append(new_hit)
    return(hits_list)


def write_to_fasta_parallel(dic_seq, num_file):
    '''Write sequences in fasta file, and separate its in several files if specified.'''
    list_seq=list(dic_seq.keys())
    sep=len(dic_seq)//num_file
    list_dic_fasta=[]
    eprint("WD " + WORKDIR)
    eprint("NF " + str(num_file))
    for num in range(num_file):
        out=open(WORKDIR + '/sgRNA'+str(num)+'.fa','w')
        list_dic_fasta.append({ 'num' : num,
                      'input_fasta' : WORKDIR + '/sgRNA' + str(num) + '.fa',
                      'results' : None
                    })
        i=0
        for seq in list_seq:
            while(i < sep):
                remove_seq=list_seq.pop()
                out.write('>' + remove_seq + '\n' + remove_seq + '\n')
                i += 1
        if num == num_file-1:
            if list_seq:
                for seq in list_seq:
                    out.write('>' + seq + '\n' + seq + '\n')
        out.close()
    return(list_dic_fasta)


def run_bowtie(organism_code,fasta_file,num):
    resultFile = WORKDIR + '/results_bowtie' + num + '.sam'
    bowtie_tab=['bowtie2','-x ' + REF_GEN_DIR + '/index2/' + organism_code + '/' + organism_code + ' -f ' + fasta_file + ' -S ' + resultFile + ' -L 13 -a --quiet ']
    subprocess.call(bowtie_tab)
    return resultFile

def add_notin_parallel(num_thread,list_fasta,organism_code,dic_seq):
    '''Launch bowtie alignments for excluded genomes and treat the results, with parallelization (ie if 4 threads are selected, then 4 bowtie will be launch at the same time, with 4 subfiles of the beginning file.
    For excluded genomes, only the sequence NOT matching with genome will be conserved. 
    '''
    def worker():
        while True:
            e = q.get()
            fasta_file=e['input_fasta']
            num_str=str(e['num'])
            resultFile = run_bowtie(organism_code,fasta_file,num_str)
            #resultFile = WORKDIR + '/results_bowtie' + num_str + '.sam'
            res=open(resultFile, 'r')
            dic_result={}
            count=0
            for l in res:
                if l[0]!='@':
                    if l.split('\t')[2]=='*':
                        count+=1
                        seq=l.split('\t')[0]
                        dic_result[seq]=dic_seq[seq]
            e['results']=dic_result
            res.close()
            q.task_done()

    q = Queue()
    for i in range(num_thread):
        t = Thread(target=worker)
        t.daemon = True
        t.start()

    for e in list_fasta:
        q.put(e)

    q.join()
    total_results={}
    for e in list_fasta:
        total_results.update(e['results'])

    return total_results


def add_in_parallel(num_thread,list_fasta,organism_code,dic_seq,genome,len_sgrna):
    '''Launch bowtie alignments for included genomes and treat the results, with parallelization (ie if 4 threads are selected, then 4 bowtie will be launch at the same time, with 4 subfiles of the beginning file.
    For included genomes, only the sequences matching exactly (no mismatch0 with genome will be conserved. 
    '''
    
    def worker():
        while True:
            e = q.get()
            fasta_file=e['input_fasta']
            num_str=str(e['num'])
            resultFile = run_bowtie(organism_code, fasta_file, num_str)
            res=open(resultFile, 'r')
            dic_result={}
            for l in res:
                if l[0]!='@':
                    if l.split('\t')[2]!='*':
                        l_split=l.split('\t')
                        cigar=l_split[5]
                        if cigar=='23M':
                            mm=l_split[-2]
                            ref=l_split[2]
                            if mm.split(':')[-1]=='23':
                                seq=l_split[0]
                                if seq not in dic_result:
                                    dic_result[seq]=dic_seq[seq]
                                if genome not in dic_result[seq]:
                                    dic_result[seq][genome]=[]
                                seq_align=l_split[9]
                                if seq != seq_align:
                                    strand='+'
                                else:
                                    strand='-'
                                start=l_split[3]
                                end=int(start)+len_sgrna-1
                                coord=ref+':'+strand+'('+start+','+str(end)+')'
                                dic_result[seq][genome].append(coord)
            e['results']=dic_result
            res.close()
            q.task_done()

    q = Queue()
    for i in range(num_thread):
        t = Thread(target=worker)
        t.daemon = True
        t.start()

    for e in list_fasta:
        q.put(e)

    q.join()
    total_results={}
    for e in list_fasta:
        total_results.update(e['results'])

    return total_results

def sort_hits(hitlist):
    '''
    Scoring of the hits found, where positive scores mean stronger sgRNA constructs.
    Complete Hit objects with score and sort the hits with this scores.
    '''
    eprint('-- Sort hits --')
    for hit in hitlist:
        for genome in hit.genomes_Dict:
            hit.score+=len(hit.genomes_Dict[genome])
    sorted_hitlist=sorted(hitlist,key=lambda hit:hit.score,reverse=True)
    return(sorted_hitlist)


def write_to_file(genomes_IN,genomes_NOT_IN,hit_list,PAM,non_PAM_motif_length):
    '''Write results in a file.
    The file is a tabulated file, with first column=sgRNA sequence, then one column for each included organisms with list of positions of the sequence in each.
    '''
    eprint('-- Write results to file --')
    responseResultFile = WORKDIR + '/results_allgenome.txt'
    output=open(responseResultFile,'w')
    not_in=True
    gi=','.join(genomes_IN)
    gni=','.join(genomes_NOT_IN)
    if gni=='':
        gni='None'
        not_in=False
    output.write('#ALL GENOMES\n#Genomes included :'+gi+' ; Genomes excluded :'+gni+'\n'+'#Parameters: PAM:'+PAM+' ; sgRNA size:'+str(non_PAM_motif_length)+'\n')
    output.write('sgRNA sequence')
    for genome_i in genomes_IN:
        output.write('\t'+genome_i)
    output.write('\n')
    for hit in hit_list:
        output.write(hit.sequence)
        for gi in genomes_IN:
            output.write('\t'+','.join(hit.genomes_Dict[gi]))
        output.write('\n')
    output.close()

def output_interface(hit_list):
    '''
    Reformat the results to create a json file. 
    It will be parsed in javascript to display it in interface.
    '''
    eprint('-- Construct results for graphical interface --')
    json_result_file=WORKDIR+'/results.json'
    #print(json_result_file)
    list_dic=[]

    for hit in hit_list:
        dic_json={'sequence':hit.sequence,'occurences':[]}
        list_coords=[]
        for genome in hit.genomes_Dict: 
            dic_coords={'org':genome,'coords':hit.genomes_Dict[genome]}
            list_coords.append(dic_coords)
        dic_json['occurences']=list_coords    
        list_dic.append(dic_json)

    
    with open(json_result_file,'w') as f: 
        json.dump(list_dic,f)


def order_for_research(list_in,list_notin,genome,dict_org_code,dist_dic,list_order):
    '''Determinate the order of research, based on the distance matrix created for all genomes. 
    The best comparisons to do are comparisons between similar genomes for the exclusion and between distant genomes for inclusion
    '''

    ref1=dict_org_code[genome][0]
    if list_in and list_notin:
        in_compare=-1
        for gi in list_in:
            ref2=dict_org_code[gi][0]
            dist=dist_dic[ref1][ref2]
            if dist > in_compare:
                in_compare=dist
                in_compare_genome=gi
        notin_compare=11
        for gni in list_notin:
            ref2=dict_org_code[gni][0]
            dist=dist_dic[ref1][ref2]
            if dist < notin_compare:
                notin_compare=dist
                notin_compare_genome=gni

        if in_compare > notin_compare :
            new_genome=in_compare_genome
            list_order.append((new_genome,'in'))
            list_in.remove(new_genome)
        else:
            new_genome=notin_compare_genome
            list_order.append((new_genome,'notin'))
            list_notin.remove(new_genome)

    elif list_in:
        in_compare=-1
        for gi in list_in:
            ref2=dict_org_code[gi][0]
            dist=dist_dic[ref1][ref2]
            if dist > in_compare:
                in_compare=dist
                in_compare_genome=gi
        new_genome=in_compare_genome
        list_order.append((new_genome,'in'))
        list_in.remove(new_genome)

    elif list_notin:
        notin_compare=11
        for gni in list_notin:
            ref2=dict_org_code[gni][0]
            dist=dist_dic[ref1][ref2]
            if dist < notin_compare:
                notin_compare=dist
                notin_compare_genome=gni
        new_genome=notin_compare_genome
        list_order.append((new_genome,'notin'))
        list_notin.remove(new_genome)
    else:
        return(list_order)

    return order_for_research(list_in,list_notin,new_genome,dict_org_code,dist_dic,list_order)

def unzip_files(list_genomes,dict_org_code):
    '''Unzip required files for reserach'''
    eprint('-- Unzip selected genomes --')
    start=time.time()
    out=open(WORKDIR+'/unzip.sh','w')
    for genome in list_genomes:
        ref=dict_org_code[genome][0] 
        out.write('tar xf '+REF_GEN_DIR+'/fasta/'+ref+'.tar.gz\ntar xf '+REF_GEN_DIR+'/index2/'+ref+'.tar.gz\n')
    out.close()
    os.system('bash ' + WORKDIR + '/unzip.sh')  
    end=time.time()
    eprint('TIME UNZIP',end-start)

def delete_used_files(list_genomes,dict_org_code):
    '''Delete unzipped files that have been used for research'''
    eprint('-- Delete selected genomes --')
    start=time.time()
    out=open(WORKDIR+'/delete.sh','w')
    for genome in list_genomes:
        ref=dict_org_code[genome][0] 
        out.write('rm -r '+REF_GEN_DIR+'/fasta/'+ref+'\n')
        out.write('rm -r '+REF_GEN_DIR+'/index2/'+ref+'\n')
    out.close()
    os.system('bash '+WORKDIR+'/delete.sh')
    end=time.time()
    eprint('TIME DELETE', end-start)   

def construction(fasta_path,PAM,non_PAM_motif_length,genomes_IN,genomes_NOT_IN,dict_org_code):

    start_time=time.time()
    start = time.time()
    num_thread=4
    num_file=4
    eprint('## Search for',len(genomes_IN),"included genomes and",len(genomes_NOT_IN),'excluded genomes with',num_thread,'thread(s) ##')

    unzip_files(genomes_IN+genomes_NOT_IN,dict_org_code)

    if len(genomes_IN)!=1:
        sorted_genomes=sort_genomes(genomes_IN,fasta_path,dict_org_code)
    else:
        sorted_genomes=genomes_IN

    if len(genomes_NOT_IN)>=1:
        sorted_genomes_notin=sort_genomes_desc(genomes_NOT_IN,fasta_path,dict_org_code)
    else:
        sorted_genomes_notin=[]

    eprint('-- RESEARCH --')    
    dic_seq=construct_in(fasta_path,sorted_genomes[0],dict_org_code[sorted_genomes[0]][0],PAM,non_PAM_motif_length)
    eprint(str(len(dic_seq))+' hits in first included genome '+sorted_genomes[0])
    list_fasta=write_to_fasta_parallel(dic_seq,num_file)
    distFile = REF_GEN_DIR + "/distance_dic.json"
    #eprint("Pickle Fiellocation is " + pickleFile)
    dist_dic = None
    with open(distFile, 'r') as f:
        dist_dic = json.load(f)
    #dist_dic=json.load(open(pickleFile, "rb" ))
    eprint('-- Determinate order for research -- ')
    list_order=order_for_research(sorted_genomes[1:],sorted_genomes_notin,sorted_genomes[0],dict_org_code,dist_dic,[])

    for i in list_order:
        genome=i[0]
        if i[1]=='notin':
            dic_seq=add_notin_parallel(num_thread,list_fasta,dict_org_code[genome][0],dic_seq)
            if len(dic_seq)==0:
                print("Program terminated&No hits remain after exclude genome "+genome)
                end_time=time.time()
                total_time=end_time-start_time
                eprint('TIME',total_time)
                delete_used_files(genomes_NOT_IN+genomes_IN,dict_org_code)
                sys.exit(1)
            eprint(str(len(dic_seq))+' hits remain after exclude genome '+genome)
            list_fasta=write_to_fasta_parallel(dic_seq,num_file)

        elif i[1]=='in':
            dic_seq=add_in_parallel(num_thread,list_fasta,dict_org_code[genome][0],dic_seq,genome,len(PAM)+non_PAM_motif_length)
            if len(dic_seq)==0:
                print("Program terminated&No hits remain after include genome "+genome)
                end_time=time.time()
                total_time=end_time-start_time
                eprint('TIME',total_time)
                delete_used_files(genomes_NOT_IN+genomes_IN,dict_org_code)
                sys.exit(1)
            eprint(str(len(dic_seq))+' hits remain after include genome '+genome)
            list_fasta=write_to_fasta_parallel(dic_seq,num_file)

    delete_used_files(genomes_IN+genomes_NOT_IN,dict_org_code)       
    print(len(dic_seq))        

    hit_list=construct_hitlist(dic_seq)

    hit_list=sort_hits(hit_list)

    ##Put results in local file for access via the interface.
    write_to_file(genomes_IN,genomes_NOT_IN,hit_list[:10000],PAM,non_PAM_motif_length)

    ##Output formatting for printing to interface
    output_interface(hit_list[:100])

def main():

    eprint("--->" + TASK_KEY)
    start_time=time.time()

    parameters = args_gestion()

    #indexs_path = parameters.rfg + '/index2'

    fasta_path  = parameters.rfg + '/fasta'
    global REF_GEN_DIR
    REF_GEN_DIR = parameters.rfg
    dict_organism_code = readOrganismCode(REF_GEN_DIR + '/genome_ref_taxid.json')  ##Keys: organism, values: genomic reference (ncbi)
   # organisms_selected,organisms_excluded,PAM,non_PAM_motif_length=args_gestion(dict_organism_code)
    organisms_selected, organisms_excluded, PAM, non_PAM_motif_length = setupApplication(parameters, dict_organism_code)
    print(','.join(organisms_excluded))
    print(TASK_KEY)
    global WORKDIR
    WORKDIR = setupWorkSpace(parameters)

    eprint('---- CSTB complete genomes ----')
    eprint('Parallelisation with distance matrix')
    construction(fasta_path,PAM,non_PAM_motif_length,organisms_selected,organisms_excluded,dict_organism_code)
    end_time=time.time()
    total_time=end_time-start_time
    eprint('TIME',total_time)
    #eprint('CUMULATIVE_LENGTH',cumulative_length(organisms_selected,dict_organism_code))
    #eprint('CUMULATIVE_LENGTH_NOT_IN',cumulative_length(organisms_excluded,dict_organism_code))

if __name__=='__main__':
    main()