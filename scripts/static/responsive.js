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

function readTextFile(file)	//Found at http://stackoverflow.com/questions/14446447/javascript-read-local-text-file
{
    var rawFile = new XMLHttpRequest();
    rawFile.open("GET", file, false);
    rawFile.onreadystatechange = function ()
    {
        if(rawFile.readyState === 4)
        {
            if(rawFile.status === 200 || rawFile.status == 0)
            {
                var allText = rawFile.responseText;
                $("#INList").html(allText)
                $("#gref").html(allText)
                $("#INList_sg").html(allText)

            }
        }
    }
    rawFile.send(null);
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

// Modif begin

function openCity(evt, cityName) {
  var i, x, tablinks;
  x = document.getElementsByClassName("city");
  for (i = 0; i < x.length; i++) {
     x[i].style.display = "none";
  }
  tablinks = document.getElementsByClassName("tablink");
  for (i = 0; i < x.length; i++) {
     tablinks[i].className = tablinks[i].className.replace(" w3-border-red", "");
  }
  document.getElementById(cityName).style.display = "block";
  evt.currentTarget.firstElementChild.className += " w3-border-red";
  $('#tree_include').jstree().show_dots();
  $('#tree_exclude').jstree().show_dots();
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
		print_seq='<td rowspan="'+line_seq+'">' + seq + '</td>'
		out=print_seq+out

	}
	var header='<tr class = "w3-light-grey"> <th> sgRNA sequence </th> <th> Organism(s) </th> <th colspan=2> Coordinates </th> </tr>'
	out=header+out
	return out

}


function displayTreeIn(suffix){
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


function selectOnTreeIn(data){
	for(m = 0, n = includeIdList.length; m < n; m++) {
		disabledExc.push(includeIdList[m].replace("j1","j2"));
	}
	$("#tree_exclude").jstree().enable_node(disabledExc);

	includeIdList = data.selected
	disabledExc = []
	for(m = 0, n = includeIdList.length; m < n; m++) {
		disabledExc.push(includeIdList[m].replace("j1","j2"));
		}
	$("#tree_exclude").jstree().disable_node(disabledExc);

	includeNameList = [];
	for(i = 0, j = includeIdList.length; i < j; i++) {
		includeNameList.push(data.instance.get_node(includeIdList[i]).text);
	};

	notSelectedIdIn=[]
	for(i=0;i<excludeIdList.length;i++){
		notSelectedIdIn.push(excludeIdList[i].replace('j2','j1'))
	}
	$('#tree_include').jstree(true).uncheck_node(notSelectedIdIn);
}

function selectOnTreeNotIn(data){
	for(m = 0, n = excludeIdList.length; m < n; m++) {
		disabledInc.push(excludeIdList[m].replace("j2","j1"));
	}
	$("#tree_include").jstree().enable_node(disabledInc);
	excludeIdList = data.selected
	disabledInc = []
	for(m = 0, n = excludeIdList.length; m < n; m++) {
		disabledInc.push(excludeIdList[m].replace("j2","j1"));
	}
	$("#tree_include").jstree().disable_node(disabledInc);
	excludeNameList = [];
	for(i = 0, j = excludeIdList.length; i < j; i++) {
		excludeNameList.push(data.instance.get_node(excludeIdList[i]).text);
	};

	notSelectedIdNotIn=[]
	for(i=0;i<includeIdList.length;i++){
		notSelectedIdNotIn.push(includeIdList[i].replace('j1','j2'))
	}
	$('#tree_exclude').jstree(true).uncheck_node(notSelectedIdNotIn);
	
}

function resetTree(){
	includeNameList = [];
	excludeNameList = [];
	$('#tree_include').jstree().close_all();
	$('#tree_include').jstree().deselect_all();
	$('#tree_exclude').jstree().close_all();
	$('#tree_exclude').jstree().deselect_all();
	$('#tree_include').jstree().search('');
	$('#tree_exclude').jstree().search('');
}

function submitTree(){
	$("#IN").hide();
	$("#NOTIN").hide();
	$("#RemoveIN").hide();
	$("#ValIN").hide();
	$("#RemoveNOTIN").hide();
	$("#ValNOTIN").hide();
	$("#submit_trees").hide();
	$("#reset_trees").hide();
	$('#ShowIN').show();
	$('#ShowNOTIN').show();
	displaySelection()
	$("#other_parameters_AllG").show();
	$("#submitbtn").show();
	
}

function displaySelection(){
	for (i=0;i<includeNameList.length;i++){
		var n=new Option(includeNameList[i]);
		$(n).html(includeNameList[i]);
		$("#InView").append(n);	//Adds the contents of 'includeNameList' array into the 'InView'
	};
	for (i=0;i<excludeNameList.length;i++){
		var n=new Option(excludeNameList[i]);
		$(n).html(excludeNameList[i]);
		$("#NotInView").append(n);
	};	
}

function resetAll(){
	$("#INList").empty();
	$("#InView").empty();
	$("#NOTINList").empty();
	$("#NotInView").empty();
	$("#other_parameters_AllG").hide();
	$("#ShowIN").hide();
	$("#NOTIN").hide();
	$("#ShowNOTIN").hide();
	$("#submitbtn").hide();
	$("#submit_trees").show();
	$("#reset_trees").show();
	$("#IN").show();
}
// Modif end

function submitSetupAllGenome(){
	GIN=JSON.stringify(includeNameList);
	GNOTIN=JSON.stringify(excludeNameList);
	MISMATCH=$("select[name='max-mismatch_AllG'] > option:selected").val();
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

$(document).ready(function(){

	readTextFile("./static/sortedgenomes.txt")	//Requires putting the sortedgenomes.txt in the 'static' folder.

	displayTreeIn('')

	displayTreeNotIn('')

	includeIdList=[]
	disabledExc=[]
	excludeIdList=[]
	disabledInc=[]
	excludeNameList=[]
	includeNameList=[]

	$('#tree_include').on("changed.jstree",function(e,data){
		selectOnTreeIn(data)
	})

	$('#tree_exclude').on("changed.jstree", function (e, data) {	// replace changed for onclicked like below
		selectOnTreeNotIn(data)

	});



//END OF TREE

	$("#rstbtn").click(resetAll)

	$("#reset_trees").click(resetTree)

	$('#search_in').val('');
	$('#search_notin').val('');

	$("#submit_trees").click(submitTree)

	$("#submitbtn").click(function(){
		$("#AllG").hide();
		$("#Waiting").show();
		submitSetupAllGenome()
		$.getJSON($SCRIPT_ROOT+'/allgenomes',{gi:GIN,gni:GNOTIN,max_mismatch:MISMATCH,pam:PAM,sgrna_length:SGRNA},
		function(data) {
			treatResults(data)
		});	

	});	
		


	//**SPECGENE CODE**//



	var includeIdList_sg = []
	var excludeIdList_sg = []
	var includeNameList_sg = []	// this is the list contain the names of IN organisms
	var excludeNameList_sg = []	// this is the list contain the names of NOT IN organisms
	var disabledExc_sg = []
	var disabledInc_sg = []

	var gref='';
	var res_specGene='';
	var Not_In_search=0;		//Holds binary information as to whether Not In search performed (1) or not (0)
	var Multiple_search=0;

	$('#changeSeq').hide()
	$('#GREF').hide()
	$('#IN_sg').hide()
	$('#NOTIN_sg').hide()
	$("#other_parameters").hide()
	$('#submit_sg').hide()
	$('#tree_sg').hide()

	$('#load-file').click(function(){
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
	});

	$('#next').click(function(){
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
		// YOLO
		$('#ShowSeq').html("<h5 class='w3-container w3-light-green'>Your query:</h5><div class='w3-margin'>"+sequence+"</div>")
		//A embellir+mettre dans le HTML?

		$('#spec_tips').hide()
		$('a[name=error-fasta]').hide()
		$("#Sequenceupload").hide()
		$('#next').hide()
		$('#list').hide()
		$('#changeSeq').show()
		$('#tree_sg').show()
	});


	$(function () {
		// search on IN tree
		var to = false;
		$('#search_in_sg').keyup(function () {
			if(to) { clearTimeout(to); }
			to = setTimeout(function () {
				var v = $('#search_in_sg').val();
				$('#tree_include_sg').jstree().search(v);
			}, 250);
		});
		$.jstree.defaults.search.show_only_matches=true;
		// tree building
		$('#tree_include_sg').jstree({
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
		$('#tree_include_sg').jstree().show_dots();
	});


	// The second tree : NOT IN
	$(function () {
		// search on NOT IN tree
		var to = false;
		$('#search_notin_sg').keyup(function () {
			if(to) { clearTimeout(to); }
			to = setTimeout(function () {
				var v = $('#search_notin').val();
				$('#tree_exclude').jstree().search(v);
			}, 250);
		});
		$.jstree.defaults.search.show_only_matches=true;
		// tree building
		$('#tree_exclude_sg').jstree({
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
		$('#tree_exclude_sg').jstree().show_dots();
	});

	$('#tree_include_sg').on("changed.jstree", function (e, data) {	// replace changed for onclicked like below
		for(m = 0, n = includeIdList_sg.length; m < n; m++) {
			disabledExc_sg.push(includeIdList_sg[m].replace("j1","j2"));
		}
		$("#tree_exclude_sg").jstree().enable_node(disabledExc_sg);
		includeIdList_sg = data.selected
		disabledExc_sg = []
		for(m = 0, n = includeIdList_sg.length; m < n; m++) {
			disabledExc_sg.push(includeIdList_sg[m].replace("j1","j2"));
		}
		$("#tree_exclude_sg").jstree().disable_node(disabledExc_sg);
		includeNameList_sg = [];
		for(i = 0, j = includeIdList_sg.length; i < j; i++) {
			includeNameList_sg.push(data.instance.get_node(includeIdList_sg[i]).text);
		};
		console.log(includeNameList_sg);
	});

	$('#tree_exclude_sg').on("changed.jstree", function (e, data) {	// replace changed for onclicked like below
		for(m = 0, n = excludeIdList_sg.length; m < n; m++) {
			disabledInc_sg.push(excludeIdList_sg[m].replace("j2","j1"));
		}
		$("#tree_include_sg").jstree().enable_node(disabledInc_sg);
		excludeIdList_sg = data.selected
		disabledInc_sg = []
		for(m = 0, n = excludeIdList_sg.length; m < n; m++) {
			disabledInc_sg.push(excludeIdList_sg[m].replace("j2","j1"));
		}
		$("#tree_include_sg").jstree().disable_node(disabledInc_sg);
		excludeNameList_sg = [];
		for(i = 0, j = excludeIdList_sg.length; i < j; i++) {
			excludeNameList_sg.push(data.instance.get_node(excludeIdList_sg[i]).text);
		};
		console.log(excludeNameList_sg);
	});

	var genomes_IN_sg=[];
	var genomes_NOTIN_sg=[]

	var out_specGene=''	//Holds the HTML Table content for output display.

	$('#submit_trees_sg').click(function(){
		genomes_IN=includeNameList_sg;
		genomes_NOTIN=excludeNameList_sg;
		$('#tree_sg').hide()
		for (i=0;i<includeNameList_sg.length;i++){
			var n=new Option(includeNameList_sg[i]);
			$(n).html(includeNameList_sg[i]);
			//$("#InView").append(n);	//Adds the contents of 'includeNameList' array into the 'InView'
		};
		for (i=0;i<excludeNameList_sg.length;i++){
			var n=new Option(excludeNameList_sg[i]);
			$(n).html(excludeNameList_sg[i]);
			//$("#NotInView").append(n);
		};
		$("#other_parameters").show();
		$("#submit_sg").show();
	});

	$('#submit_sg').click(function(){

		$("#Tabselection").hide();
		$('a[name=error-n]').hide()
		$('a[name=error-pid]').hide()
		var n_gene=$("#search-region").val()
		var percent_id=$("#percent-identity").val()
		var pam=$("select[name='pam'] > option:selected").val();
		var sgrna_length=$("select[name='sgrna-length'] > option:selected").val();
		//ERRORS
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

		$('#SpecG').hide();
		$('#Waiting').show();
		var N=JSON.stringify(n_gene);
		var PID=JSON.stringify(percent_id);
		var SEQ=JSON.stringify(final_sequence);
		var GIN=JSON.stringify(genomes_IN);
		var GNOTIN=JSON.stringify(genomes_NOTIN);
		var GREF=JSON.stringify(gref);
		var PAM=JSON.stringify(pam);
		var LENGTH=JSON.stringify(sgrna_length);

		$.getJSON($SCRIPT_ROOT+'/specific_gene',{seq:SEQ,gin:GIN,gnotin:GNOTIN,gref:GREF,n:N,pid:PID,pam:PAM,sgrna_length:LENGTH},
			function(data) {
			if (data.length==5){
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
			display_download(tag)	//Link to where the output file is held; used for downloading the results.
			//Same file used for AllG and SpecG options.
			}
			else{		//Display no matching results output.
				$("#Waiting").hide();
				$("#NoResult").show();
				infos='<p>'+data[0]+'</p> <p> '+data[1]+'</p>'
				$("#no_result").html(infos);
			}
		});
	genomes_IN=[];	//These ready the lists
	genomes_NOTIN=[];	//for the next request.
	});


});

function reloadpage() {
    location.reload();
}