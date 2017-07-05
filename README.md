# CSTB

## Requirements 

This program has been developped with Python3 (v3.6.1), and needs Flask (v0.12.1 used), BioPython (v1.69 used), Bowtie2 (v2.3.1 used) and Blastn (v2.5.0+ used) 

### Python

#### On Linux (Ubuntu) 
`sudo apt-get install python3`

#### On MacOs  
One easy solution is to use the [Homebrew package manager](http://brew.sh/). Follow the link and install homebrew via the terminal as specified. 
Then : `brew install python3`

[Python website](https://www.python.org)

It's also possible to use Python2. In this case, you have to change the global variable PYTHON_INTERPRETER to "python" instead of "python3" in webapp.py script (folder scripts). **Warning** : If you choose to use Python2, Flask and BioPython must be installed for this version.   

### Flask and BioPython 
You can use pip3 with `pip3 install flask` and `pip3 install biopython`

#### Online documentations
[BioPython](http://biopython.org/wiki/Documentation) 

[Flask](http://flask.pocoo.org/docs/0.12/)

### Blast and Bowtie2 

On Linux (Ubuntu), you can use `sudo apt-get install blast` and `sudo apt-get install bowtie2`. 
For more recent versions, check the documentations. 
#### Online documentations 

[Blast](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

### Warning 
If you choose to download requirements "manually", you have to add them to your PATH variable, so that they can be executable directly in your terminal.  

## How to launch 

### Construct a database 

The github archive does not provide genome database for local version. You can construct one with the scripts pre_treatment.py. 

For build database, you will need an ncbi assembly_summary file. This kind of files are available in NCBI FTP and starts with assembly_summary. 

**Example of files :** [All refseq genomes](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt), [All bacteria](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt) 

Once you download this type of file, launch `python3 scripts/pre_treatment.py path_to_assembly_summary_file` from CSTB folder to construct database. At the end, you must have a new folder called reference_genomes with different needed files. 



### Local version 
Go to CSTB folder and execute `bash run.command`. This script executes Flask and you can go to the local adress http://127.0.0.1:5000/ to access graphical interface. 

You can execute Flask manually with the following commands (you have to be in CSTB folder) : 

`export FLASK_APP=scripts/webapp.py`

`export FLASK_DEBUG=1` (optionnal, only useful if you plan to modify the code) 

`flask run`

And go the the local url to access graphical interface. 

## How it works 
CSTB has been created to find target sequence for CRISPr system. 
The programm has two part called All Genomes and Specific Gene. The functionment is the same, the only difference is that All Genomes search sequences in complete genomes and Specific Gene search sequences only in a given gene.  
The principle is to detect sequences that hybridises on targeted genomes and DON'T hybridise on excluded genomes (it's the novelty compared to the other programs that exist).  

### All Genomes principle 

1. Sort genomes. Included genomes are sorted by ascending size and excluded genomes by descending size. 
2. Search sequences in first included genome (the smallest) with Python regular expressions, and store the coordinates. 
3. Determines the order of search, for next comparisons. This order is based on genomes similarity. Indeed, the search is faster when genome have high similarity with an excluded genome or low similarity with an included genome. 
4. Map the current sequences against the next genome to compare with Bowtie2. If next genome is included, only sequences with exact match are conserved. If next genome is excluded, only sequences that don't map are conserved. The coordinates in the new genome are stored. Then, these sequences are mapped against the next genome and so on until there's no more genomes left. 
5. Sort and format the results.

### Specific Gene principle 

1. Search homologous of the given gene in included genomes with Blast. Store the coordinates of genes. 
2. Search sequences in given gene (with Python regular expression). 

The next steps are exactly the same as All Genomes steps 1 to 5, except that more informations are stored like if the sequence is present in the gene or not.

### Performances test 

*Coming soon...* 
