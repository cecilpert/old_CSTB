import sys
from ete3 import Tree,NCBITaxa
from functionbase import construct_dict_organism_assemblyref

"""
Last change : 5pm 11 dec. 2016

Run thois script with the commande below:
	
	python3 tax2json.py <tax_file>

The <tax_file> must be a text file witch has content like this:
	
	Buchnera aphidicola str. APS (Acyrthosiphon pisum)
	Borrelia burgdorferi B31
	Treponema denticola ATCC 35405
	<NCBI taxonomy name>
	...

The output file will be named 'json_tree.json' as json format or
flat file.

* Part of script adapted from https://gist.github.com/jhcepas/9205262
"""

def get_json(node):
	node.name = node.name.replace("'", '')    
	json = {"text": node.name}
	if node.children:
		json["children"] = []
		for ch in node.children:
			json["children"].append(get_json(ch))
	return json

if __name__ == '__main__':

	Dict=construct_dict_organism_assemblyref()	##This Dictionary contains as keys the names of all the bacterial genomes in the database.
	
	ncbi = NCBITaxa()
	allspe = []
	
	# get taxid from species name list, then save them into a list
	#file = open(sys.argv[1],'r')
	for sp in Dict.keys():
		ID = ncbi.get_name_translator([sp])
		ID = int(ID[sp][0])
		allspe.append(ID)
	#file.close()

	# Build topologic tree
	treeTopo = ncbi.get_topology(allspe)

	# Convert to Newick tree as string format
	treeNwk = treeTopo.write(format_root_node=True, format=8)

	# build tree object from Newick format string
	newTree = Tree(treeNwk, format=1)	# in = str
	
	# get nodes and leaves name
	for i in newTree.iter_descendants():
		i.name = ncbi.get_taxid_translator([int(i.name)])[int(i.name)]
	# get root's name
	newTree.name = ncbi.get_taxid_translator([int(newTree.name)])[int(newTree.name)]
	
	
	# Convert to json format
	json_tree = str(get_json(newTree)).replace("'", '"')
	
	# output to json file (or .txt format)
	output = open('./scripts/static/jsontree.json','w')
	output.write(json_tree)
	output.close()
