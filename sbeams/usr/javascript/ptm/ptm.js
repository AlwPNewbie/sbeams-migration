/*var ptm = [];
ptm[0] = 'GFIY(0.0225)IS(0.0330)A(0.0888)S(0.8474)HNPIGY(0.0083)NGIK';
ptm[1] = 'GFIY(0.0122)IS(0.0179)A(0.0364)S(0.9226)HNPIGY(0.0108)NGIK';
ptm[2] = 'S(0.9390)HNPIGY(0.0610)NGIK';
ptm[3] = 'IEEELGDKAVYAGENFHHGDKL';


function main() {
    // for testing...
    for (var pep of ptm) {
	document.body.appendChild(document.createElement("hr"));
	document.body.appendChild(document.createTextNode(pep));
	document.body.appendChild(document.createElement("br"));
	document.body.appendChild(document.createElement("br"));
	document.body.appendChild(show_prob_bars(pep));
	document.body.appendChild(document.createElement("br"));
	document.body.appendChild(document.createElement("br"));
    }
    document.body.appendChild(document.createElement("hr"));
}
*/
function show_prob_bars(ptmpep,barcolor='#ff5f00',label='ptmprob=',maxpix=30,width=10){
    var htmlpep = document.createElement("span");
    htmlpep.style.display = 'inline-grid';
    htmlpep.style.gridTemplateRows = maxpix+"px auto";
    htmlpep.style.gridAutoFlow = 'column';

    for (var str of ptmpep.split(")")) {
    	var [aas, prob] = str.split("(");

	    var bar;
	    for (var aa of aas.split("")) {
        bar = document.createElement("span");
				bar.style.display= 'inline-block';
				bar.style.height = maxpix + 'px';
				bar.style.width = width + "px";
				bar.style.background =  'transparent';
				bar.style.marginLeft =  '1px';
				htmlpep.appendChild(bar);

				var span = document.createElement("span");
				span.style.textAlign = 'center';
				span.appendChild(document.createTextNode(aa));
				htmlpep.appendChild(span);
	   }

		 if (prob) {
				var pct = (100 * Number(prob));
				bar.style.background = 'linear-gradient(to top, '+barcolor+', '+barcolor+' '+pct+'%, #f3f1e4 '+pct+'%)';
				bar.title = label+prob;
		 }
    }
    return htmlpep;
}

