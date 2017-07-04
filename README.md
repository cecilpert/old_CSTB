# CSTB

## Requirements 

This program has been developped with Python3 (v3.6.1), and needs Flask (v0.12.1 used), BioPython (v1.69 used), Bowtie2 (v2.3.1 used) and Blastn (v2.5.0+ used) 

### Python

#### On Linux (Ubuntu) 
*sudo apt-get install python3*

#### On MacOs  
one easy solution is to use the Homebrew package manager(http://brew.sh/). Follow the link and install homebrew via the terminal as specified. 
Then : *brew install python3*

##### Online tutorial 
https://wiki.python.org/moin/BeginnersGuide/Download for Python3

It's also possible to use Python2. In this case, you have to change the global variable PYTHON_INTERPRETER to "python" instead of "python3" in webapp.py script (folder scripts). **Warning** : If you choose to use Python2, Flask and BioPython must be installed for this version.   

### Flask and BioPython 
You can use pip3 with *pip3 install flask* and *pip3 install biopython*

#### Online documentations
<li>http://biopython.org/wiki/Download for Biopython </li>
<li> http://flask.pocoo.org/docs/0.12/ </li>

### Blast and Bowtie2 

On Linux (Ubuntu), you can use *sudo apt-get install blast* and *sudo apt-get install bowtie2*. 
For more recent versions, check the documentations. 
#### Online documentations 
<li> https://www.ncbi.nlm.nih.gov/books/NBK279690/ </li> 
<li> http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml </li> 

### Warning 
If you choose to download requirements "manually", you have to add them to your PATH variable, so that they can be executable directly in your terminal.  






