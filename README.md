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

Once you download this type of file, launch `python3 scripts/pre_treatment.py` from CSTB folder to construct database. At the end, you must have a new folder called reference_genomes with different needed files. 



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

*Coming soon...* 

### Specific Gene principle 

*Coming soon...* 

### Performances test 

*Coming soon...* 




