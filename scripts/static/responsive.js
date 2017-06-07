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

function display_download(file)
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
                onDownload(allText)
            }
        }
    }
    rawFile.send(null);
}

function onDownload(data) {
     var to_display=encodeURIComponent(data)
     $('#result-file').html('<a href="data:application/txt;charset=utf-8,'+to_display+'"download="results.txt">Download results</a>')                                     
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


// Modif end

$(document).ready(function(){

	readTextFile("./static/sortedgenomes.txt")	//Requires putting the sortedgenomes.txt in the 'static' folder.


    //**ALLGENOMES CODE**//

//TREE IMPLEMENTATION

	var includeIdList = []
	var excludeIdList = []		
	var includeNameList = []	// this is the list contain the names of IN organisms
	var excludeNameList = []	// this is the list contain the names of NOT IN organisms
	var disabledExc = []
	var disabledInc = []

	// The first tree : IN
	$(function () {
		// search on IN tree
		var to = false;
		$('#search_in').keyup(function () {
			if(to) { clearTimeout(to); }
			to = setTimeout(function () {
				var v = $('#search_in').val();
				$('#tree_include').jstree().search(v);
			}, 250);
		});
		$.jstree.defaults.search.show_only_matches=true;
		// tree building
		$('#tree_include').jstree({
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
		$('#tree_include').jstree().show_dots();
	});


	// The second tree : NOT IN
	$(function () {
		// search on NOT IN tree
		var to = false;
		$('#search_notin').keyup(function () {
			if(to) { clearTimeout(to); }
			to = setTimeout(function () {
				var v = $('#search_notin').val();
				$('#tree_exclude').jstree().search(v);
			}, 250);
		});
		$.jstree.defaults.search.show_only_matches=true;
		// tree building
		$('#tree_exclude').jstree({
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
		$('#tree_exclude').jstree().show_dots();
	});

	$('#tree_include').on("changed.jstree", function (e, data) {	// replace changed for onclicked like below
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
		console.log(includeNameList);
	});

	$('#tree_exclude').on("changed.jstree", function (e, data) {	// replace changed for onclicked like below
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
		console.log(excludeNameList);
	});


//END OF TREE

	//These variables are arrays
	var genomes_IN=[];
	var genomes_NOTIN=[]
	//Which store the target and excluded organisms.
	
	//This array will hold a copy of filled genomes_IN for output generation.
	var conserve_genomes_IN=[];
	
	//This array will store all the selectable organisms on initiation.
	var choices=[];
	//
	
	//This will store data callbacks from HTTP requests.
	var res='';
	var seq='';
	var resultFound=1;
	
	$("#INList option").each(function(){
	choices.push($(this).text());
	});
	
	
	$("#rstbtn").click(function(){	//Resets the interface to what it looks like when it is initially loaded.
		var genomes_IN=[];
		var genomes_NOTIN=[];
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
		for (i=0;i<choices.length;i++){
			var n=new Option(choices[i]);
			$(n).html(choices[i]);
			$("#INList").append(n);
		};
		$("#IN").show();
	});

	$("#reset_trees").click(function(){
		includeNameList = [];
		excludeNameList = [];
		$('#tree_include').jstree().close_all();
		$('#tree_include').jstree().deselect_all();
		$('#tree_exclude').jstree().close_all();
		$('#tree_exclude').jstree().deselect_all();
		$('#search_in').val('');
		$('#search_notin').val('');
		$('#tree_include').jstree().search('');
		$('#tree_exclude').jstree().search('');
	});
	$('#search_in').val('');
	$('#search_notin').val('');
	$('#seq').val('');

	
	$("#submit_trees").click(function(){
		genomes_IN=includeNameList;
		genomes_NOTIN=excludeNameList;
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
		$("#other_parameters_AllG").show();
		if(genomes_NOTIN.length==0){
			$("#MaxMisMatches").hide();
		}
		$("#submitbtn").show();
	});
	
	//Set hovertext for genomes IN selection box
	$("#INList").attr('title', 'You can select multiple options at once');
	//
	
	$("#AddIN").click(function(){
		$("#submit_trees").hide();
		$("#reset_trees").hide();
		$("#ShowIN").show();
		$("#RemoveIN").show();	//Ensures you see this after a reset
		$("#ValIN").show();		//Ensures you see this after a reset
		var names=[];
		$("#INList option:selected").each(function(){	//Takes the selected options in the list of genomes and puts them in the 'names' array
			names.push($(this).text());
		});
		for (i=0;i<names.length;i++){
			var n=new Option(names[i]);
			$(n).html(names[i]);
			$("#InView").append(n);	//Adds the contents of 'names' array into the 'InView'
			$("#INList option:contains('"+names[i]+"')").remove();
		};
	});	
	
	$("#RemoveIN").click(function(){
		var names=[];
		$("#InView option:selected").each(function(){
			names.push($(this).text());
		});
		for (i=0;i<names.length;i++){
			var n=new Option(names[i]);
			$(n).html(names[i]);
			$("#INList").append(n);
			$("#InView option:contains('"+names[i]+"')").remove();
		};
	});	
	
	$("#ValIN").click(function(){
		var foo=[];
		$("#IN").hide();
		$("#RemoveIN").hide();	//We show the list of genomes_IN but disallow remove
		$("#ValIN").hide();		//And Validate choices option buttons by hiding them.
		$("#INList option").each(function(){
			foo.push($(this).text());
		});
		$("#InView option").each(function(){
			genomes_IN.push($(this).text());
		});
		foo.sort();
		for (i=0;i<foo.length;i++){
			var n=new Option(foo[i]);
			$(n).html(foo[i]);
			$("#NOTINList").append(n);
		};
		$("#NOTIN").show();
	});
	
	$("#AddNOTIN").click(function(){
		$("#ShowNOTIN").show();
		$("#ValNOTIN").show();
		$("#RemoveNOTIN").show();
		var names=[];
		$("#NOTINList option:selected").each(function(){
			names.push($(this).text());
		});
		for (i=0;i<names.length;i++){
			var n=new Option(names[i]);
			$(n).html(names[i]);
			$("#NotInView").append(n);
			$("#NOTINList option:contains('"+names[i]+"')").remove();
		};
		});
		
	$("#NoneNOTIN").click(function(){
		$("#RemoveNOTIN").hide();
		$("#ShowNOTIN").show();
		$("#ValNOTIN").show();
	});
	
	$("#RemoveNOTIN").click(function(){
		var names=[];
		$("#NotInView option:selected").each(function(){
			names.push($(this).text());
		});
		for (i=0;i<names.length;i++){
			var n=new Option(names[i]);
			$(n).html(names[i]);
			$("#NOTINList").append(n);
			$("#NotInView option:contains('"+names[i]+"')").remove();
		};
	});	
	
	
	$("#ValNOTIN").click(function(){
		var foo=[];
		$("#NOTIN").hide();
		$("#RemoveNOTIN").hide();
		$("#ValNOTIN").hide();
		$("#NotInView option").each(function(){
			genomes_NOTIN.push($(this).text());
		});
		$("#other_parameters_AllG").show();
		if(genomes_NOTIN.length==0){
			$("#MaxMisMatches").hide();
		}
		$("#submitbtn").show();
	});
		
		
	$("#submitbtn").click(function(){
		$("#Searchoptions").hide();
		$("#AllG").hide();
		$("#allg_tips").hide();
		$("#Waiting").show();
		$("#Tabselection").hide();
		var GIN=JSON.stringify(genomes_IN);
		var GNOTIN=JSON.stringify(genomes_NOTIN);
		var MISMATCH=$("select[name='max-mismatch_AllG'] > option:selected").val();
		var PAM=$("select[name='pam_AllG'] > option:selected").val();
		var SGRNA=$("select[name='sgrna-length_AllG'] > option:selected").val();
		$.getJSON($SCRIPT_ROOT+'/allgenomes',{gi:GIN,gni:GNOTIN,max_mismatch:MISMATCH,pam:PAM,sgrna_length:SGRNA},
		function(data) {
			if (data.length==3){
				res=data[0];
				window.alert(res)
				not_in=data[1]; 
				window.alert(not_in)
				tag=data[2];
				window.alert(tag)
			}
			else {
			resultFound=0;
			}
			if(resultFound==1){
				var obj=JSON.parse(res);
				var out='<tr><th colspan=2>All hits are absent (based on Bowtie2 alignment) from excluded genomes : '+not_in+'</th></tr><tr><td colspan=2></td></th>';
				$("#Waiting").fadeOut();
				$("#Result").show();
				for (i=0;i<obj.length;i++){
					out+='<tr><th colspan=2>Construct: '+obj[i].sequence+'</th></tr>';
					out+='<tr><th>Genomes: </th><th>Coordinates: </th></tr>';
					for (j=0;j<obj[i].in.length;j++){
						out+='<tr><td>'+obj[i].in[j].org+'</td> <td>'+obj[i].in[j].coords+'</td></tr>';
					}	
				}	
			$("#ResTable").html(out);
			display_download("./static/results"+tag+".txt")	//Link to where the output file is held; used for downloading the results.
			//Same file used for AllG and SpecG options.
			}
			else{		//Display no matching results output.
				$("#Waiting").hide();
				$("#Result").show();
				$("#output").text(data[0]+'\n'+data[1]);
			}
		});
		genomes_IN=[];	//These ready the lists
		genomes_NOTIN=[];	//for the next request.
	});


	//**SPECGENE CODE**// 

	var sequence='';
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

	$("#INList_sg option").each(function(){
		choices.push($(this).text());
	});
	
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
		$('#GREF').show()
	});


	$('#valGREF').click(function(){
		gref=$('#gref').val()
		if (gref==''){
			alert('You have to select an origin genome')
			return
		}
		else if(gref.length>1){
			alert("You can't select more than one origin genome")
			return
		}
		else{
			$('#GREF').hide()
			$('#IN_sg').show()
		}	
		document.getElementById("showGref").innerHTML = '<b> Your origin genome : </b>'+gref;
		$("#INList_sg option:contains('"+gref+"')").remove()
	});

	$("#AddIN_sg").click(function(){
		$("#ShowIN_sg").show();
		$("#RemoveIN_sg").show();	//Ensures you see this after a reset
		$("#ValIN_sg").show();		//Ensures you see this after a reset
		var names=[];
		$("#INList_sg option:selected").each(function(){
			names.push($(this).text());
		});
		for (i=0;i<names.length;i++){
			var n=new Option(names[i]);
			$(n).html(names[i]);
			$("#InView_sg").append(n);
			$("#INList_sg option:contains('"+names[i]+"')").remove();
		};
	});

	$("#NoneIN_sg").click(function(){
		$("#RemoveIN_sg").hide();
		$("#ShowIN_sg").show();
		$("#ValIN_sg").show();
	});

	$("#RemoveIN_sg").click(function(){
		var names=[];
		$("#InView_sg option:selected").each(function(){
			names.push($(this).text());
		});
		for (i=0;i<names.length;i++){
			var n=new Option(names[i]);
			$(n).html(names[i]);
			$("#INList_sg").append(n);
			$("#InView_sg option:contains('"+names[i]+"')").remove();
		};
	});

	$("#ValIN_sg").click(function(){
		var foo=[];
		$("#IN_sg").hide();
		$("#RemoveIN_sg").hide();	//We show the list of genomes_IN but disallow remove
		$("#ValIN_sg").hide();		//And Validate choices option buttons by hiding them.
		$("#INList_sg option").each(function(){
			foo.push($(this).text());
		});
		$("#InView_sg option").each(function(){
			genomes_IN.push($(this).text());
		});
		foo.sort();
		for (i=0;i<foo.length;i++){
			var n=new Option(foo[i]);
			$(n).html(foo[i]);
			$("#NOTINList_sg").append(n);
		};
		$("#NOTIN_sg").show();
	});	

	$("#AddNOTIN_sg").click(function(){
		$("#ShowNOTIN_sg").show();
		$("#ValNOTIN_sg").show();
		$("#RemoveNOTIN_sg").show();
		var names=[];
		$("#NOTINList_sg option:selected").each(function(){
			names.push($(this).text());
		});
		for (i=0;i<names.length;i++){
			var n=new Option(names[i]);
			$(n).html(names[i]);
			$("#NotInView_sg").append(n);
			$("#NOTINList_sg option:contains('"+names[i]+"')").remove();
		};
	});
	
	$("#NoneNOTIN_sg").click(function(){
		$("#RemoveNOTIN_sg").hide();
		$("#ShowNOTIN_sg").show();
		$("#ValNOTIN_sg").show();
	});
	
	$("#RemoveNOTIN_sg").click(function(){
		var names=[];
		$("#NotInView_sg option:selected").each(function(){
			names.push($(this).text());
		});
		for (i=0;i<names.length;i++){
			var n=new Option(names[i]);
			$(n).html(names[i]);
			$("#NOTINList_sg").append(n);
		$("#NotInView_sg option:contains('"+names[i]+"')").remove();
		};
	});

	$("#ValNOTIN_sg").click(function(){
		var foo=[];
		$("#NOTIN_sg").hide();
		$("#RemoveNOTIN_sg").hide();
		$("#ValNOTIN_sg").hide();
		$("#NotInView_sg option").each(function(){
		genomes_NOTIN.push($(this).text());
		});
		$("#n").show()
		$('#other_parameters').show()
		$("#submit_sg").show();
	});


	var out_specGene=''	//Holds the HTML Table content for output display.
	

	$('#submit_sg').click(function(){

		$("#Tabselection").hide();
		$('a[name=error-n]').hide()
		$('a[name=error-pid]').hide()
		var n_gene=$("#search-region").val()
		var percent_id=$("#percent-identity").val()
		var max_mismatch=$("select[name='max-mismatch'] > option:selected").val();
		var pam=$("select[name='pam'] > option:selected").val();
		var sgrna_length=$("select[name='sgrna-length'] > option:selected").val();
		var mismatch_og=$("select[name='max-mismatch-og'] > option:selected").val();
		var mismatch_nin=$("select[name='max-mismatch-notin'] > option:selected").val();
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
		var MISMATCH=JSON.stringify(max_mismatch);
		var PAM=JSON.stringify(pam);
		var LENGTH=JSON.stringify(sgrna_length);
		 
		$.getJSON($SCRIPT_ROOT+'/specific_gene',{seq:SEQ,gin:GIN,gnotin:GNOTIN,gref:GREF,n:N,pid:PID,max_mismatch:MISMATCH,pam:PAM,sgrna_length:LENGTH,max_mm_og:mismatch_og,max_mm_notin:mismatch_nin},function(data){
			if (data.length==3){	//If the Spec Gene search program terminated, the return data type is not 'string' but 'object'.
				res_specGene=data[0]
				tag=data[1]
				var obj_specGene=JSON.parse(res_specGene);
				$("#Waiting").fadeOut();
				$("#Result").show();

				if (obj_specGene[0].otherorgs!="No search"){
					Multiple_search=1
				}

				if (obj_specGene[0].not_in!="No search"){	//=if a Not in search has been requested
					Not_In_search=1
				}

				if (Multiple_search==1){ //there is more than 1 genome 

					for (i=0;i<obj_specGene.length;i++){
						var other_orgs=obj_specGene[i].otherorgs.split("%")
						out_specGene+='<tr><th>Construct: '+obj_specGene[i].refsequence+'</th></tr>';
						out_specGene+='<tr><th>Genomes: '+'</th><th>Coordinates: '+'</th></tr>';
						out_specGene+='<tr><td>'+obj_specGene[i].reforg+'</td><td>'+obj_specGene[i].on_off+'</td><tr>';

						
						
						for (j=0;j<other_orgs.length;j++){
							focus_org=other_orgs[j]
							var occs=obj_specGene[i].other_on_off[0][focus_org]	//Brackets [] allows access to name stored in variable focus_org.
							//SPLIT HERE TO ACCOMODATE '%' (multiple occurences)
							out_specGene+='<tr><td>'+focus_org+'</td><td>'+occs+'</td><tr>';
						}

						if (Not_In_search==1){	//Not In Results Display
							out_specGene+='<tr><th>Not In</th></tr>';
							var Not_Ins=obj_specGene[i].not_in_genomes.split("%")
							var Not_in_write=''
							for(k=0;k<Not_Ins.length;k++){
								if (obj_specGene[i].not_in[0][Not_Ins[k]]=="Absent"){
									Not_in_write='<b>Absent</b>'
								}
								else{
									Not_in_write='<b>Present</b>, '+obj_specGene[i].not_in[0][Not_Ins[k]]
								}
								out_specGene+='<tr><td>'+Not_Ins[k]+'</td><td>'+Not_in_write+'</td><tr>'
							}
						}
					}
				}

				else{ //there is only one genome 

					for (i=0;i<obj_specGene.length;i++){
						out_specGene+='<tr><th>Construct: '+obj_specGene[i].refsequence+'</th></tr>';
						out_specGene+='<tr><td>Reference: '+obj_specGene[i].reforg+'\t Coordinates: '+obj_specGene[i].on_off+'</td></tr>';
						
						if (Not_In_search==1){	//Not In Results Display
							out_specGene+='<tr><th>Not In</th></tr>';
							var Not_Ins=obj_specGene[i].not_in_genomes.split("%")
							var Not_in_write=''
							for(k=0;k<Not_Ins.length;k++){
								
								if (obj_specGene[i].not_in[0][Not_Ins[k]]==""){
									Not_in_write='<b>Present, exact match</b>'
								}
								else if (obj_specGene[i].not_in[0][Not_Ins[k]]=="No matches"){
									Not_in_write='<b>Absent</b>'
								}
								else{
									Not_in_write='<b>Present</b>, with substitutions '+obj_specGene[i].not_in[0][Not_Ins[k]]
								}
								out_specGene+='<tr><td>'+Not_Ins[k]+"\t"+Not_in_write+'</td></tr>'
							}
						}
					}
				
				}
				$("#ResTable").html(out_specGene);
				display_download("./static/results"+tag+".txt")	//Link to where the output file is held; used for downloading the results.
				//Same file used for AllG and SpecG options.

			}
			else{		//Display no matching results output.
				$("#Waiting").hide();
				$("#Result").show();
				$("#output").text(data[0]+'\n'+data[1]);
			}
		});	
	genomes_IN=[];	//These ready the lists
	genomes_NOTIN=[];	//for the next request.
	});
		
	
});	

function reloadpage() {
    location.reload();
}