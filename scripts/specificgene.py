from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio import motifs
from Bio.Alphabet import IUPAC
from functionbase import reverse_complement,construct_dict_organism_assemblyref,build_expression,reverse_complement_mismatch,timestamp,cumulative_length,corrected_mismatch
import argparse, os, sys, time, re, subprocess, uuid, json
from allgenomes import readOrganismCode, setupWorkSpace, sort_genomes, sort_genomes_desc, find_sgRNA_seq, order_for_research
from queue import Queue
from threading import Thread

TASK_KEY = str(uuid.uuid1())
WORKDIR = None
REF_GEN_DIR = None


class Hit():
    def __init__(self,sequence,in_dic):
        self.sequence=sequence
        self.genomes_Dict=in_dic
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
    parser.add_argument('-pam',help='PAM motif for sgRNA', required=True)
    parser.add_argument('-sl',help='sgRNA length (without PAM motif)', required=True, type=int)
    parser.add_argument("-rfg",metavar="<str>",help="The path to the reference genome folder")
    parser.add_argument("-cah",metavar="<str>",help="The path to cache folder for the webservice")
    args = parser.parse_args()

    return args 

def setupApplication(parameters,dict_organism_code): 
    eprint(dict_organism_code)
    organisms_selected=parameters.gi.split('+')
    eprint(organisms_selected)
    organisms_excluded=parameters.gni.split('+')
    organisms_selected=[i for i in organisms_selected if i in dict_organism_code]
    organisms_excluded=[i for i in organisms_excluded if i in dict_organism_code]
    eprint(organisms_selected)
    non_PAM_motif_length=int(parameters.sl)

    return parameters.seq,organisms_selected,organisms_excluded,non_PAM_motif_length,parameters.n,parameters.ip,parameters.pam


def eprint(*args, **kwargs):
    '''For printing only on the terminal'''
    print(*args, file=sys.stderr, **kwargs)


### METHODS FOR THE RESEARCH ###

def blast_to_find_all_genes(query_seq, genomes, identity_percentage_min,dic_genome):
    '''Launch blastn between query sequence and reference genome and all in genomes (separately).
    For origin genome, will only save the genes position, not the subject sequence (the query sequence is conserved for origin genome)
    For other genomes IN, will take the subject sequence found by blast as sequence for the organism and the sgRNAs will be searched on it. It will also return gene's genomic position
    The hits must have a %id greater than the %id given by user 
    If there is several blast hits, we only take the first (as defined by Blast)
    '''

    dic_genes={} 
    size_gene=len(query_seq)
    file_query=WORKDIR+'/blast_request'
    with open(file_query,'w') as f:
        f.write(query_seq)
    f.close()    
    for i in range(len(genomes)): 
        eprint("Searching for sequence in",genomes[i])
        ref=dic_genome[genomes[i]][0]
        db_path = REF_GEN_DIR + "/fasta/" + ref + '/' + ref + "_genomic.fna" #peut-etre mieux faire sur les cds parce que avec le génome complet pas d'infos sur nom gene, description etc...
        blast_command = "blastn -db " + db_path + " -query " + file_query + " -outfmt 5"
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
                                
                        dic_genes[genomes[i]]=(gene_start,gene_end)

                    else: #other in genomes, now we take the hit subject as seq for sgRNA searching  
                        if first_hsp.sbjct_start > first_hsp.sbjct_end: #reverse strand 
                            gene_start=first_hsp.sbjct_end-1 #blast is on 1-based and our program need 0-based 
                            gene_end=first_hsp.sbjct_start-1
                            strand='-'
                        else:
                            gene_start=first_hsp.sbjct_start-1
                            gene_end=first_hsp.sbjct_end-1  
                            strand='+'
                        dic_genes[genomes[i]]=(gene_start,gene_end)                                                 
    return dic_genes 

                                           
def sgRNA_on_gene(start_gene, end_gene, start_sgRNA, end_sgRNA):
    '''Check if a sgRNA is in a gene (with both genomics positions)'''
    on_gene = False
    if start_gene <= start_sgRNA and end_gene >= end_sgRNA:
        on_gene = True
    return on_gene


def score_and_sort(hitlist):
    '''Score and sort the results of an analysis.
    For now, the score for 1 sgRNA is the sum of all mismatchs + 1 if it's present in an excluded organism (no matter if it's exact match or not)  
    Fill the attribute score of Hit objects.
    So best sgRNA is the one that have lowest score. 
    ''' 
    number_on_gene=0
    for hit in hitlist:
        dic_on_gene={}
        eprint(hit.sequence,hit.genomes_Dict)
        count_on_gene=0
        for genome in hit.genomes_Dict: 
            dic_on_gene[genome]=0
            count_genomes_on_gene=0
            list_coords=hit.genomes_Dict[genome]
            for coord in list_coords: 
                if coord.split(':')[-1]=='OnGene': 
                    count_on_gene+=1
                    dic_on_gene[genome]+=1
        hit.score=count_on_gene  
        count=0
        for genome in dic_on_gene: 
            if dic_on_gene[genome]>=1: 
                count+=1
        if count==len(hit.genomes_Dict): 
            number_on_gene+=1
    
    print(number_on_gene)            
    sorted_hitlist=sorted(hitlist,key=lambda hit:hit.score,reverse=True)    
    return sorted_hitlist
  

def find_sgRNA_given_sequence(fasta_path, organism, organism_code, PAM, non_PAM_motif_length):
    '''
    Same function as construct_seq_dict except the dictionnary will not be like value=list of coordinates. It add the information about organism with dictionnary values like : dictionnary with key=organism and value=list of coordinates in this organism
    '''
    fasta_file=fasta_path + '/' + organism_code + '/' + organism_code +'_genomic.fna'

    sgRNA=''
    for i in range(non_PAM_motif_length):
        sgRNA+='N'
    sgRNA+=PAM

    seq_dict={}

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

def search_sgRNA_in_query_seq(query_seq,non_PAM_motif_length, PAM): 

    sgRNA=''
    for i in range(non_PAM_motif_length):
        sgRNA+='N'
    sgRNA+=PAM    

    dic_seq={}

    seq_list_forward=find_sgRNA_seq(query_seq,reverse_complement(sgRNA))

    for indice in seq_list_forward:

        end=indice+len(PAM)+non_PAM_motif_length
        seq=reverse_complement(query_seq[indice:end])
        dic_seq[seq]={}

    return dic_seq    

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

def add_in_parallel(num_thread,list_fasta,organism_code,dic_seq,genome,len_sgrna,gene_coords):
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
                                dic_result[seq][genome]=[]
                                seq_align=l_split[9]
                                if seq != seq_align:
                                    strand='+'
                                else:
                                    strand='-'
                                start=l_split[3]
                                end=int(start)+len_sgrna-1
                                if sgRNA_on_gene(gene_coords[0],gene_coords[1],int(start),int(end)): 
                                    info_gene='OnGene'
                                else:     
                                    info_gene='OffGene'
                                coord=ref+':'+strand+'('+start+','+str(end)+'):'+info_gene
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

def setupWorkSpace(parameters):
    '''Create the work folder where all temporary or results files will be stored'''
    workFolder = parameters.cah + '/' + TASK_KEY
    os.mkdir(workFolder)
    return workFolder   

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

def do_search(query_seq, n, genome_list, dict_org_code, not_in_genome_list,
              identity_percentage_min,pam,sgrna_length):
    '''Launch the research with all parameters given by user. Principal function, it will call other functions to do the research'''

    num_threads=2
    fasta_path=REF_GEN_DIR+'/fasta'
    len_all_sgrna=sgrna_length+len(pam)
    distFile = REF_GEN_DIR + "/distance_dic.json"
    with open(distFile, 'r') as f:
        dist_dic = json.load(f)
    f.close()


    unzip_files(genome_list+not_in_genome_list,dict_org_code)

    dic_genes=blast_to_find_all_genes(query_seq, genome_list, identity_percentage_min,dict_org_code)

    dic_seq=search_sgRNA_in_query_seq(query_seq,sgrna_length,pam)

    sorted_genomes_in=sort_genomes(genome_list,fasta_path,dict_org_code)

    if not_in_genome_list: 
        sorted_genomes_notin=sort_genomes_desc(not_in_genome_list,fasta_path,dict_org_code)
    else: 
        sorted_genomes_notin=[]    

    first_genome=sorted_genomes_in[0]
    list_fasta=write_to_fasta_parallel(dic_seq,num_threads)

    dic_seq=add_in_parallel(num_threads,list_fasta,dict_org_code[first_genome][0],dic_seq,first_genome,len_all_sgrna,dic_genes[first_genome])
    list_fasta=write_to_fasta_parallel(dic_seq,num_threads)

    eprint(dic_seq)
    list_order=order_for_research(sorted_genomes_in[1:],sorted_genomes_notin,sorted_genomes_in[0],dict_org_code,dist_dic,[])

    for i in list_order:
        genome=i[0]
        if i[1]=='notin':
            dic_seq=add_notin_parallel(num_threads,list_fasta,dict_org_code[genome][0],dic_seq)
            if len(dic_seq)==0:
                print("Program terminated&No hits remain after exclude genome "+genome)
                end_time=time.time()
                total_time=end_time-start_time
                eprint('TIME',total_time)
                delete_used_files(not_in_genome_list+genome_list,dict_org_code)
                sys.exit(1)
            eprint(str(len(dic_seq))+' hits remain after exclude genome '+genome)
            list_fasta=write_to_fasta_parallel(dic_seq,num_threads)

        elif i[1]=='in':
            dic_seq=add_in_parallel(num_threads,list_fasta,dict_org_code[genome][0],dic_seq,genome,len_all_sgrna,dic_genes[genome])
            if len(dic_seq)==0:
                print("Program terminated&No hits remain after include genome "+genome)
                end_time=time.time()
                total_time=end_time-start_time
                eprint('TIME',total_time)
                delete_used_files(not_in_genome_list+genome_list,dict_org_code)
                sys.exit(1)
            eprint(str(len(dic_seq))+' hits remain after include genome '+genome)
            list_fasta=write_to_fasta_parallel(dic_seq,num_threads)

    delete_used_files(not_in_genome_list+genome_list,dict_org_code)       
    print(len(dic_seq))   


    hit_list=construct_hitlist(dic_seq)

    hit_list=score_and_sort(hit_list)

    write_to_file(genome_list,not_in_genome_list,hit_list[:10000],pam,sgrna_length)

    ##Output formatting for printing to interface
    output_interface(hit_list[:100]) 
   

def main():
    start_time=time.time()

    parameters = args_gestion()

    global REF_GEN_DIR
    REF_GEN_DIR = parameters.rfg

    dict_organism_code = readOrganismCode(REF_GEN_DIR + '/genome_ref_taxid.json')

    global WORKDIR
    WORKDIR = setupWorkSpace(parameters)

    query_seq,organisms_selected,organisms_excluded,non_PAM_motif_length,n,identity_percentage_min,pam=setupApplication(parameters,dict_organism_code)
    print(','.join(organisms_excluded))
    print(TASK_KEY)
    do_search(query_seq, n, organisms_selected, dict_organism_code, organisms_excluded, identity_percentage_min, pam, non_PAM_motif_length)
    
    end_time=time.time()
    total_time=end_time-start_time
    eprint('TIME',total_time)
    #eprint('CUMULATIVE_LENGTH',cumulative_length(genome_list,dic_genome))

if __name__ == "__main__":
    main()
