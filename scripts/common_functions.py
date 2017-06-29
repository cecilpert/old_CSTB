import os,sys,json,re,subprocess
from Bio import SeqIO
from queue import Queue
from threading import Thread

class Hit():
    def __init__(self,sequence,in_dic):
        self.sequence=sequence
        self.genomes_Dict=in_dic
        self.score=0


def eprint(*args, **kwargs):
    '''For printing only on the terminal'''
    print(*args, file=sys.stderr, **kwargs) 


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


def find_sgRNA_seq(seq,pam):
    '''
    Uses Regular expression matching of the PAM motif to the reference genome to get the start
    positions (0-based) of each match
    '''
    list_seq=[]
    reg_exp = build_expression(pam)
    indices = [m.start() for m in re.finditer('(?=' + reg_exp + ')', seq, re.I)]
    return indices

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

def unzip_files(list_genomes,dict_org_code):
    '''Unzip required files for reserach'''
    eprint('-- Unzip selected genomes --')
    out=open(WORKDIR+'/unzip.sh','w')
    for genome in list_genomes:
        ref=dict_org_code[genome][0] 
        out.write('tar xf '+REF_GEN_DIR+'/fasta/'+ref+'.tar.gz\ntar xf '+REF_GEN_DIR+'/index2/'+ref+'.tar.gz\n')
    out.close()
    os.system('bash ' + WORKDIR + '/unzip.sh')  

def write_to_fasta_parallel(dic_seq, num_file):
    '''Write sequences in fasta file, and separate its in several files if specified.'''
    list_seq=list(dic_seq.keys())
    sep=len(dic_seq)//num_file
    list_dic_fasta=[]
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

def treat_bowtie_notin(resultFile,e,dic_seq):    
    res=open(resultFile, 'r')
    dic_result={}
    for l in res:
        if l[0]!='@':
            if l.split('\t')[2]=='*':
                seq=l.split('\t')[0]
                dic_result[seq]=dic_seq[seq]
    e['results']=dic_result
    res.close()
   

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
            treat_bowtie_notin(resultFile,e,dic_seq)
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
    out=open(WORKDIR+'/delete.sh','w')
    for genome in list_genomes:
        ref=dict_org_code[genome][0] 
        out.write('rm -r '+REF_GEN_DIR+'/fasta/'+ref+'\n')
        out.write('rm -r '+REF_GEN_DIR+'/index2/'+ref+'\n')
    out.close()
    os.system('bash '+WORKDIR+'/delete.sh')


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
                           
def readJsonDic(path):
    '''Load the dictionnary that contains genomes name associated with their references'''
    with open(path) as json_data:
        d = json.load(json_data)
    return d   

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
