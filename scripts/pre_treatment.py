import os, sys
from Bio import SeqIO
from ete3 import NCBITaxa
import pickle 
from functionbase import reverse_complement
import re, json
PROXY_SETTINGS = None
import argparse



def main():
    global PROXY_SETTINGS
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="List of bacteteria in \"assembly-like\" format")  
    parser.add_argument("--proxy", help="proxy settings")
    

    args = parser.parse_args()
    if args.proxy:
        PROXY_SETTINGS = args.proxy
        os.environ["http_proxy"] = args.proxy
        os.environ["https_proxy"] = args.proxy
      
    #ref_bacteria='more_genomes/assembly_summary_bacteria_500.txt'
    #dic_taxid=dic_download(ref_bacteria)
   # print ("Using " + sys.argv[1] + " as input")
    
    if not os.path.exists("./reference_genomes"):
        os.makedirs("./reference_genomes")
        os.makedirs("./reference_genomes/fasta")
        os.makedirs("./reference_genomes/index2")
        os.makedirs("./reference_genomes/pre_calculate")

        
    else :
        print('A folder names "reference_genomes" alrady exist please move it')
        sys.exit(1)
    
    dic_taxid=dic_download(args.input)
    eliminate_plasmides(dic_taxid)
#test(dic_taxid)
    index_bowtie_blast(dic_taxid)
    pre_calculate(dic_taxid)
    distance_matrix(dic_taxid)
    
    #json_tree()

class Lineage: 
    def __init__(self):
        self.species="No specie"
        self.genus="No genus"   
        self.family="No family"
        self.order="No order"
        self.classe='No class'
        self.phylum="No phylum"

def dic_download(ref_bacteria):
    print('DOWNLOAD')
    out=open('to_download.sh','w')
    taxfile=open('scripts/taxfile.txt','w')
    count=0
    f=open(ref_bacteria,'r')
    cmd=''
    dic_taxid={}
    dic_ref={}
    for l in f: 
        count+=1
        l_split=l.rstrip().split('\t')
        name=l_split[7]
        short_ref=l_split[0]
        taxfile.write(name+'\n')
        ftp_link=l_split[19]
        ref=ftp_link.split('/')[-1]
        taxid=l_split[5]
        dic_taxid[ref]=taxid
        dic_ref[name+' '+short_ref]=[ref,taxid]
        ftp_suffix=ftp_link.split('/')[-1]+'_genomic.fna'
        genome_link='https://'+ftp_link.split('//')[1]+'/'+ftp_suffix
        if ftp_suffix not in os.listdir('reference_genomes/fasta'): 
            cmd += 'curl '
            if PROXY_SETTINGS:
                cmd += '--proxy ' + PROXY_SETTINGS + ' '
            cmd+='--remote-name --remote-time '+genome_link+'.gz\n'
            cmd+='mv '+ftp_suffix+'.gz reference_genomes/fasta/\n'
            cmd+='gunzip reference_genomes/fasta/'+ftp_suffix+'.gz\n'
                
    out.write(cmd)
    out.close()     
    taxfile.close()
    f.close()
    json.dump(dic_ref,open('reference_genomes/genome_ref_taxid.json','w'),indent=4)
    os.system('bash to_download.sh')
    os.system('rm to_download.sh')
    return dic_taxid 
        

def eliminate_plasmides(list_ref):     
    print('PLASMIDES')
    count=0
    for organism_code in list_ref: 
        count+=1
        print(count)
        fasta_file='reference_genomes/fasta/' + organism_code +'_genomic.fna'
        genome_seqrecord=[]
        for seq_record in SeqIO.parse(fasta_file,'fasta'):
            if 'plasmid' not in seq_record.description: 
                genome_seqrecord.append(seq_record)
        SeqIO.write(genome_seqrecord,'reference_genomes/fasta/'+organism_code+'_genomic.fna','fasta')

def index_bowtie_blast(list_ref): 
    print('INDEX')
    out=open('index.sh','w')
    for ref in list_ref: 
        if ref+'.1.bt2' not in os.listdir('reference_genomes/index2'): 
            cmd='bowtie2-build reference_genomes/fasta/'+ref+'_genomic.fna reference_genomes/index2/'+ref
            out.write(cmd+'\n')
        if ref+'_genomic.fna.nin' not in os.listdir('reference_genomes/fasta'): 
            cmd='makeblastdb -in reference_genomes/fasta/'+ref+'_genomic.fna -dbtype nucl'  
            out.write(cmd+'\n')     
    out.close()   
    os.system('bash index.sh')
    os.system('rm index.sh')  

def store_positions(seq_list_forward,seq_list_reverse,organism_code,PAM,non_PAM_motif_length):
    '''
    Same function as construct_seq_dict except the dictionnary will not be like value=list of coordinates. It add the information about organism with dictionnary values like : dictionnary with key=organism and value=list of coordinates in this organism
    '''
    fasta_file='reference_genomes/fasta/' + organism_code +'_genomic.fna'
    genome_seqrecord=next(SeqIO.parse(fasta_file, 'fasta'))
    genome_seq=genome_seqrecord.seq
    seq_dict={}
    for indice in seq_list_forward: 
        end=indice+len(PAM)+non_PAM_motif_length
        seq=genome_seq[indice:end].reverse_complement()
        seq=str(seq)
        if seq not in seq_dict: 
            seq_dict[seq]=[]
        seq_dict[seq].append('+('+str(indice+1)+','+str(end)+')')    

    for indice in seq_list_reverse: 
        end=indice+len(PAM)+non_PAM_motif_length
        seq=genome_seq[indice:end]
        seq=str(seq)
        if seq not in seq_dict: 
            seq_dict[seq]=[]
        seq_dict[seq].append('-('+str(indice+1)+','+str(end)+')')
    
    pickle.dump(seq_dict,open("reference_genomes/pre_calculate/"+organism_code+"_dicpos.pic","wb"))


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

def pre_calculate(list_ref): 
    print('PRECALCULATE')
    count=0
    for ref in list_ref: 
        count+=1
        print(count)
        if ref+'_dicpos.pic' not in os.listdir('reference_genomes/pre_calculate'): 
            list_forward,list_reverse=find_sgRNA(ref,'NGG',20)
            store_positions(list_forward,list_reverse,ref,'NGG',20)

def create_lineage_objects(dic_tax):
    ncbi=NCBITaxa()
    dic_lineage={}
    count=0
    for ref in dic_tax: 
        lineage_object=Lineage()
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
    return dic_lineage  

def distance_dic(dic_lineage): 
    dic={}
    for ref1 in dic_lineage: 
        dic[ref1]={}
        for ref2 in dic_lineage: 
            if ref1==ref2: 
                dic[ref1][ref2]=0
            elif dic_lineage[ref1][0].species==dic_lineage[ref2][0].species: 
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

def distance_matrix(dic_taxid): 
    print('DISTANCE MATRIX')
    dic_lineage=create_lineage_objects(dic_taxid)
    print("s2")
    dist_dic=distance_dic(dic_lineage)
    print("s3")
    pickle.dump(dist_dic, open( "reference_genomes/distance_dic.pickle", "wb" ) )

def test(dic_taxid): 
    print('CHECK FOR :')
    for ref in dic_taxid: 
        count=0
        for seqrecord in SeqIO.parse('reference_genomes/fasta/'+ref+'_genomic.fna','fasta'): 
            count+=1    
        if count>1: 
            print(ref)

def json_tree(): 
    print('JSON TREE')
    os.system('python3 scripts/tax2json.py')            




if __name__ == "__main__":
    main()
