function sendQuery(div,url,args) {	
	
	var params="";
	var scroll=0;
	if(url.indexOf('?')<0){params='?';}
	var arglist = args.split(',');
	for (var i=0; i<arglist.length;i++){	
		if(arglist[i] == 'scroll')scroll=1;
		var field=document.getElementById(arglist[i]);	
		if(!field){continue;}
		if(typeof(field[0])=='object'){
			var optlist = new Array();
			for(var k=0;k<field.options.length;k++){
				if(field.options[k].selected){
					optlist.push(field.options[k].value);
				}
			}
			params += "&" + arglist[i] + "=" + optlist.join('|');
		}
		else if(field.type == 'textarea'){
			params += "&" + arglist[i] + "=" + field.value.replace('\n','::::');	
		}
		else if(field.type == 'radio'){
			if( field.checked ){
				params += "&" + field.name + "=" + field.value;	
			}
		}
		else{
			params += "&" + arglist[i] + "=" + document.getElementById(arglist[i]).value;
		}
	}
		
	
	var resultDiv = document.getElementById(div);
	resultDiv.innerHTML='';
	url+=params;
	//alert(params);
	new Ajax.Request(url, {
		method:'get',
		onLoading: function() {
//			$(resultDiv).update('<img style="margin:-10em 0em;" src="/images/loading.gif"/>');
//			if (window.XMLHttpRequest){
//				resultDiv.innerHTML = '<img style="margin:-10em 0em;" src="/images/loading.gif"/>'; 
//			}
		},
		onComplete: function(transport) {
			$(resultDiv).update(transport.responseText);
		},
		onSuccess: function(transport) {
			$(resultDiv).update(transport.responseText);
		}
	});
	if(scroll){
		window.location.hash=div;
	}
}






function loadXMLDoc(div,url)
{
	if (window.XMLHttpRequest)
	  {// code for IE7+, Firefox, Chrome, Opera, Safari
	  xmlhttp=new XMLHttpRequest();
	  xmlhttp.open("GET",url,false);
	  xmlhttp.send(null);
	  }
	else
	  {// code for IE6, IE5
	  xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
	  xmlhttp.open("GET",url,false);
	  xmlhttp.send();
	  }
	document.getElementById(div).innerHTML=xmlhttp.responseText;
}
function gbrowse(){
	var organism=document.getElementById('organism').value;
	var contig=document.getElementById('contig');
	window.location="/gb2/gbrowse/pipeline?genome="+organism+"&name=" + contig.value;
}
