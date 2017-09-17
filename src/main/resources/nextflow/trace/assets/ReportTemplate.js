// JavaScript used to power the Nextflow Report Template output.
window.data_byprocess = {};
$(function() {
  // Add in some extra fields and collect stats by process
  for(i=0; i<window.data['trace'].length; i++){
    var myString = "something format_abc";
    var name_regex = /([^\(]+) \((.*)\)$/g;
    var rematch = name_regex.exec(window.data['trace'][i]['name']);
    if(rematch){
      window.data['trace'][i]['process_name'] = rematch[1];
      window.data['trace'][i]['tag'] = rematch[2];
      if(!window.data_byprocess.hasOwnProperty(rematch[1])){
        window.data_byprocess[rematch[1]] = [];
      }
      window.data_byprocess[rematch[1]].push(window.data['trace'][i]);
    } else {
      window.data['trace'][i]['process_name'] = window.data['trace'][i]['name'];
      window.data['trace'][i]['tag'] = '-';
    }
  }

  // Plot histograms of resource usage
  var cpu_data = [];
  var mem_data = [];
  var time_data = [];
  for(var pname in window.data_byprocess){
    if (window.data_byprocess.hasOwnProperty(pname)) {
      var c = [];
      var m = [];
      var t = [];
      // TODO: Hardcoded to only handle Gb
      for (var i = 0; i < window.data_byprocess[pname].length; i ++) {
        c[i] = parseInt(window.data_byprocess[pname][i]['pct_cpu']);
        m[i] = parseInt(window.data_byprocess[pname][i]['vmem']);
        t[i] = parseInt(window.data_byprocess[pname][i]['duration']);
      }
      cpu_data.push({y: c, name: pname, type:'box', boxpoints: 'all', jitter: 0.3});
      mem_data.push({y: m, name: pname, type:'box', boxpoints: 'all', jitter: 0.3});
      time_data.push({y: t, name: pname, type:'box', boxpoints: 'all', jitter: 0.3});
    }
  }
  Plotly.newPlot('cpuplot', cpu_data, { title: 'Workflow CPU Usage', yaxis: {title: '% single core CPU usage'} });
  Plotly.newPlot('memplot', mem_data, { title: 'Workflow Memory Usage', yaxis: {title: 'Memory (bytes)'} });
  Plotly.newPlot('timeplot', time_data, { title: 'Workflow Task Durations', yaxis: {title: 'Duration (milliseconds)'} });

  // Build the trace table
  for(i=0; i<window.data['trace'].length; i++){
    var status = window.data['trace'][i]['status'];
    if(status == 'COMPLETED'){ status = '<span class="badge badge-success">COMPLETED</span>'; }
    if(status == 'CACHED'){ status = '<span class="badge badge-secondary">CACHED</span>'; }
    if(status == 'FAILED'){ status = '<span class="badge badge-danger">FAILED</span>'; }
    $('#tasks_table tbody').append('<tr>'+
      '<td>' + window.data['trace'][i]['process_name'] + '</td>'+
      '<td>' + window.data['trace'][i]['tag'] + '</td>'+
      '<td>' + status + '</td>'+
      '<td>' + window.data['trace'][i]['hash'] + '</td>'+
      '<td>' + window.data['trace'][i]['name'] + '</td>'+
      '<td>' + window.data['trace'][i]['task_id'] + '</td>'+
      '<td>' + window.data['trace'][i]['realtime'] + '</td>'+
      '<td>' + window.data['trace'][i]['pct_cpu'] + '</td>'+
      '<td>' + window.data['trace'][i]['submit'] + '</td>'+
      '<td>' + window.data['trace'][i]['vmem'] + '</td>'+
      '<td>' + window.data['trace'][i]['native_id'] + '</td>'+
      '<td>' + window.data['trace'][i]['exit'] + '</td>'+
      '<td>' + window.data['trace'][i]['duration'] + '</td>'+
      '<td>' + window.data['trace'][i]['wchar'] + '</td>'+
      '<td>' + window.data['trace'][i]['rchar'] + '</td>'+
      '<td>' + window.data['trace'][i]['rss'] + '</td>'+
    +'</tr>');
  }
  $('#tasks_table').DataTable({
    "lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]]
  });

});