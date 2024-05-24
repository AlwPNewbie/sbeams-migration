/*$.getScript("../../usr/javascript/plotly/plotly-2.27.0.min.js", function() {
   alert("Script loaded but not necessarily executed.");
});*/


function draw_twoLines() {
	var trace = {
		 x: x,
		 y: y,
		 mode: 'lines+markers',
		 name: 'Sample Data'
	};
  var data=[tract];
	var layout = { title: title };
	var plotDiv = document.getElementById('plotDiv');
	Plotly.newPlot(plotDiv, data, layout);
}

function draw_oneBar() {
  var trace = {
     x: x,
     y: y,
     type: 'bar'
  };
  var data=[trace];
  var layout = { title: title };
  var plotDiv = document.getElementById('plotDiv');
  Plotly.newPlot(plotDiv, data, layout);
}

function draw_histogram() {
  var data = [];
  for (i=0; i<x.length; i++){
    var name = '';
    if (trace_names.length >= i+1){
      name = trace_names[i];
    }
		var trace = {
			 x: x[i],
			 type: 'histogram',
       bingroup: '1',
       name: name
		};
    data.push(trace);
  };
  var plotDiv = document.getElementById('plotDiv');
  Plotly.newPlot(plotDiv, data, layout);
}

function draw_scatterPlot(){
  var trace = {
    x: x,
    y: y,
    name: title,
    //opacity: 0.5,
    type: 'scatter',
    mode: 'markers',
    marker: { size: 4 },
    text: data_label
  };
  var data = [trace];
  var plotDiv = document.getElementById('plotDiv');
  Plotly.newPlot(plotDiv, data, layout);

	/*myPlot.on('plotly_click', function( data ){
		 for(var i=0; i < data.points.length; i++){
			 var idx = data.points[i]['pointNumber'];
			 var seq = data.points[i]['data']['text'][idx];
			 var url = mouseover_url;
			 url = url.replace('$regex_tag',seq);
			 window.open(url);
		 };
	});
   */
}

function draw_heatmap(){
  var trace={
		//z: [[1, 20, 30], [20, 1, 60], [30, 60, 1]],
    // x: ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'],
    //y: ['Morning', 'Afternoon', 'Evening'],
    z:z,
    x:x,
    y:y,
		type: 'heatmap'
	};
  var data = [trace];
  var plotDiv = document.getElementById('plotDiv');
	Plotly.newPlot('plotDiv', data, layout);
}
