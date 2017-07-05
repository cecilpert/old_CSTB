// last modif at 16 dec 2016 01:18
function verifyFasta(seq){
	var error='no error'
	var authorized=['A','T','C','G','a','t','c','g']
	var nbre_seq=0
	var seq_finale=''
	var sequence_split=seq.split('\n')
	for(i=0;i<sequence_split.length;i++){
		if (sequence_split[i][0]=='>'){
			nbre_seq+=1
			if(nbre_seq>1){
				error='More than 1 sequence'
				return [error,0]
			}
		}
		else{
			for (j=0;j<sequence_split[i].length;j++){
				if (!(authorized.includes(sequence_split[i][j]))){
					error="Wrong sequence. Only nucleotide characters authorized"
					return [error,0]
				}
				else{
					seq_finale+=sequence_split[i][j]
				}
			}

		}
	}
	return [error,seq_finale]
}

function loadFile(id){
			if ( ! window.FileReader ) {
				return alert( 'FileReader API is not supported by your browser.' );
			}
			var $i = $(id), // Put file input ID here
				input = $i[0];
			if ( input.files && input.files[0] ) {
				file = input.files[0]; // The file
				fr = new FileReader(); // FileReader instance

				fr.readAsText( file );
				fr.onload = function () {
					// Do stuff on onload, use fr.result for contents of file
					sequence=fr.result
					$('#seq').val(sequence)
				};
			}
			else {
				// Handle errors here
				alert( "File not selected or browser incompatible." )
			}
}

function display_download(tag) {
    onDownload(tag)
    return;
}

function onDownload(data) {
    //var to_display=encodeURIComponent(data)
    //$('#result-file').html('<a href="data:application/txt;charset=utf-8,'+to_display+'"download="results.txt">Download results</a>')
    $('#result-file').html('<a href="download/' + data +'" >Download results</a>')
}

function writeResults(obj){
	var out=''
	for (i=0;i<obj.length;i++){
		line_seq=0
		seq=obj[i].sequence
		number_genomes=obj[i].occurences.length
		for(j=0;j<number_genomes;j++){
			line_genome=0
			org=obj[i].occurences[j].org
			org_split=org.split(' ')
			ref_org=org_split.pop()
			number_ref=obj[i].occurences[j].all_ref.length
			for(k=0;k<number_ref;k++){
				ref=obj[i].occurences[j].all_ref[k].ref
				number_coords=obj[i].occurences[j].all_ref[k].coords.length
				line_genome+=number_coords
				line_seq+=number_coords
				for(l=0;l<number_coords;l++){
					coord=obj[i].occurences[j].all_ref[k].coords[l]
					print_coord='<td>'+coord+'</td></tr>'
					out=print_coord+out
				}
				print_ref='<td rowspan="'+number_coords+'">' + ref + '</td>'
				out=print_ref+out
			}
			print_org='<td rowspan="'+line_genome+'">' + org_split.join(' ') + '</td>'
			out=print_org+out
		}
		print_seq='<td rowspan="'+line_seq+'" valign="top">' + seq + '</td>'
		out=print_seq+out

	}
	var header='<tr class = "w3-light-grey"> <th> sgRNA sequence </th> <th> Organism(s) </th> <th colspan=2> Coordinates </th> </tr>'
	out=header+out
	return out

}


function displayTreeIn(suffix){
	includeIdList=[]
	disabledExc=[]

	var to = false;
	$('#search_in'+suffix).keyup(function () {
		if(to) { clearTimeout(to); }
		to = setTimeout(function () {
			var v = $('#search_in'+suffix).val();
			$('#tree_include'+suffix).jstree().search(v);
		}, 250);
	});
	$.jstree.defaults.search.show_only_matches=true;
	// tree building

	$('#tree_include'+suffix).jstree({
		"core" : {
			'data' : { "url" : "./static/jsontree.json", "dataType" : "json"},
			'themes' : [
				{"dots" : true}
			],
			'animation' : false
		},
		"checkbox" : {
			"keep_selected_style" : false
		},
		"plugins" : [ "wholerow" , "checkbox" , "search" , "dnd" , "types"]
	});
	$('#tree_include'+suffix).jstree().show_dots();
}

function displayTreeNotIn(suffix){
	excludeIdList=[]
	disabledInc=[]
	var to = false;
	$('#search_notin'+suffix).keyup(function () {
		if(to) { clearTimeout(to); }
		to = setTimeout(function () {
			var v = $('#search_notin'+suffix).val();
			$('#tree_exclude'+suffix).jstree().search(v);
		}, 250);
	});
	$.jstree.defaults.search.show_only_matches=true;
	// tree building

	$('#tree_exclude'+suffix).jstree({
		"core" : {
			'data' : { "url" : "./static/jsontree.json", "dataType" : "json"},
			'themes' : [
				{"dots" : true}
			],
			'animation' : false
		},
		"checkbox" : {
			"keep_selected_style" : false
		},
		"plugins" : [ "wholerow" , "checkbox" , "search" , "dnd" , "types"]
	});
	$('#tree_exclude'+suffix).jstree().show_dots();
}


function selectOnTreeIn(data,suffix){
	if(suffix=='_sg'){
		a='j3'
		b='j4'
	}
	else{
		a='j1'
		b='j2'
	}
	
	for(m = 0, n = includeIdList.length; m < n; m++) {
		disabledExc.push(includeIdList[m].replace(a,b));
	}
	$("#tree_exclude"+suffix).jstree().enable_node(disabledExc);

	includeIdList = data.selected
	disabledExc = []
	for(m = 0, n = includeIdList.length; m < n; m++) {
		disabledExc.push(includeIdList[m].replace(a,b));
		}
	$("#tree_exclude"+suffix).jstree().disable_node(disabledExc);

	includeNameList = [];
	for(i = 0, j = includeIdList.length; i < j; i++) {
		includeNameList.push(data.instance.get_node(includeIdList[i]).text);
	};

	notSelectedIdIn=[]
	for(i=0;i<excludeIdList.length;i++){
		notSelectedIdIn.push(excludeIdList[i].replace(b,a))
	}
	$('#tree_include'+suffix).jstree(true).uncheck_node(notSelectedIdIn);
}


function selectOnTreeNotIn(data,suffix){
	if(suffix=='_sg'){
		a='j3'
		b='j4'
	}
	else{
		a='j1'
		b='j2'
	}

	for(m = 0, n = excludeIdList.length; m < n; m++) {
		disabledInc.push(excludeIdList[m].replace(b,a));
	}
	$("#tree_include"+suffix).jstree().enable_node(disabledInc);
	excludeIdList = data.selected
	disabledInc = []
	for(m = 0, n = excludeIdList.length; m < n; m++) {
		disabledInc.push(excludeIdList[m].replace(b,a));
	}
	$("#tree_include"+suffix).jstree().disable_node(disabledInc);
	excludeNameList = [];
	for(i = 0, j = excludeIdList.length; i < j; i++) {
		excludeNameList.push(data.instance.get_node(excludeIdList[i]).text);
	};

	notSelectedIdNotIn=[]
	for(i=0;i<includeIdList.length;i++){
		notSelectedIdNotIn.push(includeIdList[i].replace(a,b))
	}
	$('#tree_exclude'+suffix).jstree(true).uncheck_node(notSelectedIdNotIn);
	
}

function resetTree(suffix){
	$('#tree_include'+suffix).jstree().close_all();
	$('#tree_include'+suffix).jstree().deselect_all();
	$('#tree_exclude'+suffix).jstree().close_all();
	$('#tree_exclude'+suffix).jstree().deselect_all();
	$('#tree_include'+suffix).jstree().search('');
	$('#tree_exclude'+suffix).jstree().search('');
	$('#search_in'+suffix).val('')
	$('#search_notin'+suffix).val('')
	includeIdList=[]
	disabledExc=[]
	excludeIdList=[]
	disabledInc=[]
	excludeNameList=[]
	includeNameList=[]
}

function submitTree(){
	$('#tree_ag').hide()
	$("#submit_trees").hide();
	$("#reset_trees").hide();
	$('#list_selection').show();
	$('#ShowIN').show()
	$('#ShowNOTIN').show()
	displaySelection('')
}

function submitTreeSG(){
	//$('#ShowSeq').hide()
	//$('#tree').hide()
	$('#submit_trees_sg').hide()
	$('#reset_trees_sg').hide()
	$('#ShowIN_sg').show()
	$('#ShowNOTIN_sg').show()
	$('#list_selection_sg').show()
	displaySelection('_sg')
}

function displaySelection(suffix){
	for (i=0;i<includeNameList.length;i++){
		node=includeIdList[i]
		if ($('#tree_include'+suffix).jstree().is_leaf(node)){
			n=new Option(includeNameList[i]);
			$(n).html(includeNameList[i]);
			$("#InView"+suffix).append(n);	//Adds the contents of 'includeNameList' array into the 'InView'
		}
		
	};
	for (i=0;i<excludeNameList.length;i++){
		node=excludeIdList[i]
		if ($('#tree_exclude'+suffix).jstree().is_leaf(node)){
			n=new Option(excludeNameList[i]);
			$(n).html(excludeNameList[i]);
			$("#NotInView"+suffix).append(n);	//Adds the contents of 'includeNameList' array into the 'InView'
		}
	};	
}

function confirmSelection(suffix){
	$('#confirm_y'+suffix).hide()
	$('#confirm_n'+suffix).hide()
	$('#other_parameters'+suffix).show()
	$('#submitbtn'+suffix).show()
}

function resetSelection(suffix){
	$('#list_selection'+suffix).hide()
	$("#ShowIN"+suffix).hide();
	$("#ShowNOTIN"+suffix).hide();
	$('#tree_ag').show();
	$("#submit_trees"+suffix).show();
	$("#reset_trees"+suffix).show();
	clearListView(suffix)
}

function clearListView(suffix){
	length_in=document.getElementById('InView'+suffix).options.length
	length_notin=document.getElementById('NotInView'+suffix).options.length
	for(i=0;i<length_in;i++){
		document.getElementById('InView'+suffix).options[0].remove()
	}
	for(j=0;j<length_notin;j++){
		document.getElementById('NotInView'+suffix).options[0].remove()
	}
}
// Modif end

function submitSetupAllGenome(){
	$('#Tabselection').hide()
	$('#allg_tips').hide()
	$('#list_selection').hide()
	$('#other_parameters').hide()

	$('#Waiting').show()
	GIN=JSON.stringify(includeNameList);
	GNOTIN=JSON.stringify(excludeNameList);
	PAM=$("select[name='pam_AllG'] > option:selected").val();
	SGRNA=$("select[name='sgrna-length_AllG'] > option:selected").val();

}

function treatResults(data){

	$("#Waiting").hide();
	

	if (data.length==4){

		$('#Result').show()
		res=data[0];
		not_in=data[1];
		tag=data[2];
		number_hits=data[3]

		var obj=JSON.parse(res);

		var infos='<p>' +number_hits + ' hits have been found for this research.' ;
		if (parseInt(number_hits)>100){
			infos+='. Only the best 100 are written below. Download the result file to get more hits.'
		} 
		if (parseInt(number_hits)>10000){
			infos+=' (only the best 10 000 are written to this file).</p>'
		}
		if (not_in!=''){
			infos+='<p> All hits are absent (based on Bowtie2 alignment) from excluded genome(s) : '+not_in;
		}
		else{
			infos+='<p> No excluded genomes selected.</p>'
		}
		out=writeResults(obj)
		$('#infos').html(infos)	
		$("#ResTable").html(out);
		display_download(tag)
	}	
	else{
		$("#NoResult").show();
		infos='<p>'+data[0]+'</p> <p> '+data[1]+'</p>'
		$("#no_result").html(infos);
	}

}


function treatFastaFile(){
	$('a[name=error-load]').hide()
	fastaname=$("#fasta-file").val()

	if(fastaname==''){
		$('a[name=error-load]').show()
		$('a[name=error-load]').html('No file selected')
		return
	}
	else{
		loadFile('#fasta-file')
	}	
}

function displaySequence(){
	sequence=$('#seq').val()
	error_fasta=verifyFasta(sequence)
	if(error_fasta[0]!="no error"){
		$('a[name=error-fasta]').html(error_fasta[0])
		return
	}
	final_sequence=error_fasta[1]
	if(final_sequence==''){
		$('a[name=error-fasta]').html('Empty sequence')
		return
	}
	$('#ShowSeq').html("<h5 class='w3-container w3-light-green'>Your query:</h5><div class='w3-margin'>"+sequence+"</div>")
	$('#spec_tips').hide()
	$('a[name=error-fasta]').hide()
	$("#Sequenceupload").hide()
	$('#next').hide()
	$('#list').hide()
	$('#changeSeq').show()
	$('#tree_sg').show()

}

function submitSpecificGene(){
	$("#Tabselection").hide();
	$('a[name=error-n]').hide()
	$('a[name=error-pid]').hide()
	n_gene=$("#search-region").val()
	percent_id=$("#percent-identity").val()
	pam=$("select[name='pam'] > option:selected").val();
	sgrna_length=$("select[name='sgrna-length'] > option:selected").val();	
	error_gestion()
	$('#spec_tips').hide();
	$('#ShowSeq').hide();
	$('#tree').hide()
	$('#list_selection_sg').hide();
	$('#other_parameters_sg').hide();
	$('#Waiting').show();
	var N=JSON.stringify(n_gene);
	var PID=JSON.stringify(percent_id);
	var SEQ=JSON.stringify(final_sequence);
	var GIN=JSON.stringify(includeNameList);
	var GNOTIN=JSON.stringify(excludeNameList);
	var PAM=JSON.stringify(pam);
	var LENGTH=JSON.stringify(sgrna_length);
	$.getJSON($SCRIPT_ROOT+'/specific_gene',{seq:SEQ,gin:GIN,gnotin:GNOTIN,n:N,pid:PID,pam:PAM,sgrna_length:LENGTH},
		function(data) {
			treatResultsSG(data)
	})
}

function treatResultsSG(data){
	if (data.length==5){
		resultFound=1
		res=data[0];
		not_in=data[1];
		tag=data[2];
		number_hits=data[3]
		number_on_gene=data[4]
	}
	else {
		resultFound=0;
	}
	if(resultFound==1){
		var obj=JSON.parse(res);
		var infos='<p>' + number_hits + ' hits have been found for this research. ' ;
		if (parseInt(number_hits)>100){
			infos+='Only the best 100 are written below. Download the result file to get more hits. '
		} 
		if (parseInt(number_hits)>10000){
			infos+='(only the best 10 000 are written to this file). '
		}
		infos+=number_on_gene+' of this hits hybridises at least one time with the subject gene (or an homologous). </p>'
		if (not_in!=''){
			infos+='<p>All hits are absent (based on Bowtie2 alignment) from excluded genome(s) : '+not_in+'.</p>';
		}
		else{
			infos+='<p>No excluded genomes selected. </p>'
		}
		out=writeResults(obj)
		$("#Waiting").fadeOut();
		$("#Result").show();
	
		$('#infos').html(infos)
		$("#ResTable").html(out);
		display_download(tag)	
	}
	else{		//Display no matching results output.
		$("#Waiting").hide();
		$("#NoResult").show();
		infos='<p>'+data[0]+'</p> <p> '+data[1]+'</p>'
		$("#no_result").html(infos);
	}
}

function error_gestion(){
	var errors=false
	try{
		if (n_gene=='')throw "format error";
		else if(n_gene<0)throw "can't be negative";
		else if(n_gene>0 && n_gene<=parseInt(sgrna_length)+parseInt(pam.length))throw "can't be smaller than sgRNA length"
		else if(n_gene>=final_sequence.length)throw "can't be larger than sequence length"
	}
	catch(err){
		$('a[name=error-n]').show()
		$('a[name=error-n]').html(err)
		errors=true
	}

	try{
		if (percent_id=='')throw "format error";
		else if(percent_id<0 || percent_id>100)throw"must be between 0 and 100";
	}
	catch(err){
		$('a[name=error-pid]').show()
		$('a[name=error-pid]').html(err)
		errors=true
	}
	if (errors==true){
		return
	}
}

function setupAll(){
	$('#AG_click').show()
	$('#allg_tips').hide()
	$('#tree_ag').hide()
	$('#list_selection').hide()
	$('#AG_click2').hide()
	$('#SG_click2').hide()
	$('#footer').hide()
	$('#spec_tips').hide()
	$('#list_selection_sg').hide()
	$('#Sequenceupload').hide()
	$('#tree').hide()
	$('#next').hide()
	$('#other_parameters_sg').hide()
	$('#ShowSeq').hide()
	displayTreeIn('')
	displayTreeNotIn('')
	displayTreeIn('_sg')
	displayTreeNotIn('_sg')
}

function setupAllGenome(){
	$('#footer').show()
	$('#allg_tips').show()
	$('#tree_ag').show()
	$('#SG_click').hide()
	$('#SG_click2').show()
	$('#AG_click').show()
	$('#AG_click2').hide()
	$('#search_in').val('')
	$('#search_notin').val('')
	$('#spec_tips').hide()

}

function setupSpecificGene(){
	$('#footer').show()
	$('#allg_tips').hide()
	$('#tree_ag').hide()
	$('#AG_click').hide()
	$('#AG_click2').show()
	$('#SG_click2').hide()
	$('#SG_click').show()
	$('#spec_tips').show()
	$('#Sequenceupload').show()
	$('#seq').val('')
	$('#next').show()

}


$(document).ready(function(){

	setupAll()

	$('#AG_click').click(function(){
		setupAll()
		setupAllGenome()
		resetTree('')
	})

	$('#AG_click2').click(function(){
		setupAll()
		setupAllGenome()
		resetTree('')
	})

	$('#tree_include').on("changed.jstree",function(e,data){
		selectOnTreeIn(data,'')
	})

	$('#tree_exclude').on("changed.jstree", function (e, data) {	// replace changed for onclicked like below
		selectOnTreeNotIn(data,'')
	})

	$('#reset_trees').click(function(){
		resetTree('')
	})

	$('#submit_trees').click(submitTree)

	$('#confirm_y').click(function(){
		confirmSelection('')
	})

	$('#confirm_n').click(function(){
		resetSelection('')
	})

	$('#submitbtn').click(function(){
		submitSetupAllGenome()
		$.getJSON($SCRIPT_ROOT+'/allgenomes',{gi:GIN,gni:GNOTIN,pam:PAM,sgrna_length:SGRNA},
		function(data) {
			treatResults(data)
		})

	})



	$('#SG_click').click(function(){
		setupAll()
		setupSpecificGene()
	})


	$('#SG_click2').click(function(){
		setupAll()
		setupSpecificGene()
	})

	$('#load-file').click(treatFastaFile)

	$('#next').click(function(){
		displaySequence()
		$('#ShowSeq').show()
		$('#tree').show()
		resetTree('_sg')
	}) 

	$('#tree_include_sg').on("changed.jstree",function(e,data){
		selectOnTreeIn(data,'_sg')
	})
	$('#tree_exclude_sg').on("changed.jstree",function(e,data){
		selectOnTreeNotIn(data,'_sg')
	})


	$('#reset_trees_sg').click(function(){
		resetTree('_sg')
	})

	$('#submit_trees_sg').click(submitTreeSG)

	$('#confirm_y_sg').click(function(){
		confirmSelection('_sg')
	})

	$('#confirm_n_sg').click(function(){
		resetSelection('_sg')
	})


	$('#submitbtn_sg').click(function(){
		submitSpecificGene()
	})




})

function reloadpage() {
    location.reload();
}
	

