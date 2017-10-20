// JavaScript used to power the Nextflow Report Template output.
window.data_byprocess = {};
$(function() {

  // Script block clicked
  $('#tasks_table').on('click', '.script_block', function(e){
    e.preventDefault();
    $(this).toggleClass('short');
  });

  // Completed date from now
  var completed_date = moment( $('#workflow_complete').text(), "ddd MMM DD HH:mm:ss .* YYYY" );
  if(completed_date.isValid()){
    $('#completed_fromnow').html('completed ' + completed_date.fromNow() + ', ');
  }

  // Collect stats by process
  for(i=0; i<window.data['trace'].length; i++){
    var proc = window.data['trace'][i]['process']
    if(!window.data_byprocess.hasOwnProperty(proc)){
      window.data_byprocess[proc] = [];
    }
    window.data_byprocess[proc].push(window.data['trace'][i]);
  }

  // Plot histograms of resource usage
  var cpu_data = [];
  var pct_cpu_used_data = [];
  var mem_data = [];
  var pct_mem_used_data = [];
  var time_data = [];
  var pct_time_used_data = [];
  var readwrite_data = [];
  var readwrite_hasdata = false;
  for(var pname in window.data_byprocess){
    if (window.data_byprocess.hasOwnProperty(pname)) {
      var c = [];
      var pc = [];
      var m = [];
      var pm = [];
      var t = [];
      var pt = [];
      var d = [];
      for (var i = 0; i < window.data_byprocess[pname].length; i ++) {
        c[i] = parseInt(window.data_byprocess[pname][i]['%cpu']);
        pc[i] = (parseInt(window.data_byprocess[pname][i]['%cpu']) / (parseInt(window.data_byprocess[pname][i]['cpus']) * 100)) * 100;
        m[i] = parseInt(window.data_byprocess[pname][i]['vmem']) / 1000000000;
        pm[i] = (parseInt(window.data_byprocess[pname][i]['vmem']) / parseInt(window.data_byprocess[pname][i]['memory'])) * 100
        t[i] = moment.duration( parseInt(window.data_byprocess[pname][i]['duration']) ).asMinutes();
        pt[i] =  (parseInt(window.data_byprocess[pname][i]['duration']) / parseInt(window.data_byprocess[pname][i]['time'])) * 100
        d[i] = (parseInt(window.data_byprocess[pname][i]['read_bytes']) + parseInt(window.data_byprocess[pname][i]['write_bytes'])) / 1000000000;
        if (d[i] > 0){ readwrite_hasdata = true; }
      }
      cpu_data.push({y: c, name: pname, type:'box', boxpoints: 'all', jitter: 0.3});
      pct_cpu_used_data.push({y: pc, name: pname, type:'box', boxpoints: 'all', jitter: 0.3});
      mem_data.push({y: m, name: pname, type:'box', boxpoints: 'all', jitter: 0.3});
      pct_mem_used_data.push({y: pm, name: pname, type:'box', boxpoints: 'all', jitter: 0.3});
      time_data.push({y: t, name: pname, type:'box', boxpoints: 'all', jitter: 0.3});
      pct_time_used_data.push({y: pt, name: pname, type:'box', boxpoints: 'all', jitter: 0.3});
      readwrite_data.push({y: d, name: pname, type:'box', boxpoints: 'all', jitter: 0.3});
    }
  }
  Plotly.newPlot('pctcpuplot', pct_cpu_used_data, { title: '% Requested CPU Used', yaxis: {title: '% Allocated CPUs Used', range: [0, 100]} });
  Plotly.newPlot('pctmemplot', pct_mem_used_data, { title: '% Requested Memory Used', yaxis: {title: '% Allocated Memory Used', range: [0, 100]} });
  Plotly.newPlot('pcttimeplot', pct_time_used_data, { title: '% Requested Time Used', yaxis: {title: '% Allocated Time Used', range: [0, 100]} });
  if(readwrite_hasdata){
    Plotly.newPlot('readwriteplot', readwrite_data, { title: 'Disk Read+Write', yaxis: {title: 'Read bytes + write (gb)'} });
  } else {
    $('#readwriteplot_div').hide();
  }
  // Only plot tabbed plots when shown
  $('#cpuplot_tablink').on('shown.bs.tab', function (e) {
    if($('#cpuplot').is(':empty')){
      Plotly.newPlot('cpuplot', cpu_data, { title: 'CPU Usage', yaxis: {title: '% single core CPU usage'} });
    }
  });
  $('#memplotlot_tablink').on('shown.bs.tab', function (e) {
     if($('#memplot').is(':empty')){
       Plotly.newPlot('memplot', mem_data, { title: 'Memory Usage', yaxis: {title: 'Memory (gb)'} });
     }
  });
  $('#timeplot_tablink').on('shown.bs.tab', function (e) {
     if($('#timeplot').is(':empty')){
       Plotly.newPlot('timeplot', time_data, { title: 'Task Durations', yaxis: {title: 'Duration (minutes)'} });
     }
  });

  // Build the trace table
  function make_duration(ms){
    if (ms == '-' || ms == 0){
      return ms;
    }
    return moment.duration( parseInt(ms) ).asMinutes().toFixed(2);
  }
  function make_date(ms){
    if (ms == '-' || ms == 0){
      return ms;
    }
    return moment( parseInt(ms) ).format();
  }
  for(i=0; i<window.data['trace'].length; i++){
    // Use a nice label for the status
    var status = window.data['trace'][i]['status'];
    if(status == 'COMPLETED'){ status = '<span class="badge badge-success">COMPLETED</span>'; }
    if(status == 'CACHED'){ status = '<span class="badge badge-secondary">CACHED</span>'; }
    if(status == 'ABORTED'){ status = '<span class="badge badge-danger">ABORTED</span>'; }
    if(status == 'FAILED'){ status = '<span class="badge badge-danger">FAILED</span>'; }
    // Trim whitespace properly from the script
    var script = '';
    var lines = window.data['trace'][i]['script'].split("\n");
    var ws_re = /^(\s+)/g;
    var flws_match = ws_re.exec(lines[1]);
    if(flws_match == null){
      script = window.data['trace'][i]['script'];
    } else {
      for(var j=0; j<lines.length; j++){
        script += lines[j].replace(new RegExp('^'+flws_match[1]), '').replace(/\s+$/,'') + "\n";
      }
    }
    $('#tasks_table tbody').append('<tr>'+
      '<td>' + window.data['trace'][i]['task_id'] + '</td>'+
      '<td>' + window.data['trace'][i]['process'] + '</td>'+
      '<td>' + window.data['trace'][i]['tag'] + '</td>'+
      '<td>' + status + '</td>'+
      '<td><code>' + window.data['trace'][i]['hash'] + '</code></td>'+
      '<td>' + window.data['trace'][i]['cpus'] + '</td>'+
      '<td>' + window.data['trace'][i]['%cpu'] + '</td>'+
      '<td>' + window.data['trace'][i]['memory'] + '</td>'+
      '<td>' + window.data['trace'][i]['%mem'] + '</td>'+
      '<td>' + window.data['trace'][i]['vmem'] + '</td>'+
      '<td>' + window.data['trace'][i]['rss'] + '</td>'+
      '<td>' + window.data['trace'][i]['peak_vmem'] + '</td>'+
      '<td>' + window.data['trace'][i]['peak_rss'] + '</td>'+
      '<td>' + make_duration(window.data['trace'][i]['time']) + '</td>'+
      '<td>' + make_duration(window.data['trace'][i]['duration']) + '</td>'+
      '<td>' + make_duration(window.data['trace'][i]['realtime']) + '</td>'+
      '<td><pre class="script_block short"><code>' + script.trim() + '</code></pre></td>'+
      '<td>' + window.data['trace'][i]['exit'] + '</td>'+
      '<td>' + make_date(window.data['trace'][i]['submit']) + '</td>'+
      '<td>' + make_date(window.data['trace'][i]['start']) + '</td>'+
      '<td>' + make_date(window.data['trace'][i]['complete']) + '</td>'+
      '<td>' + window.data['trace'][i]['rchar'] + '</td>'+
      '<td>' + window.data['trace'][i]['wchar'] + '</td>'+
      '<td>' + window.data['trace'][i]['syscr'] + '</td>'+
      '<td>' + window.data['trace'][i]['syscw'] + '</td>'+
      '<td>' + window.data['trace'][i]['read_bytes'] + '</td>'+
      '<td>' + window.data['trace'][i]['write_bytes'] + '</td>'+
      '<td>' + window.data['trace'][i]['native_id'] + '</td>'+
      '<td>' + window.data['trace'][i]['name'] + '</td>'+
      '<td>' + window.data['trace'][i]['module'] + '</td>'+
      '<td>' + window.data['trace'][i]['container'] + '</td>'+
      '<td>' + window.data['trace'][i]['disk'] + '</td>'+
      '<td>' + window.data['trace'][i]['attempt'] + '</td>'+
      '<td><samp>' + window.data['trace'][i]['scratch'] + '</samp></td>'+
      '<td><samp>' + window.data['trace'][i]['workdir'] + '</samp></td>'+
    +'</tr>');
  }
  $('#tasks_table').removeAttr('width').DataTable({
    "lengthMenu": [[25, 50, 100, -1], [25, 50, 100, "All"]],
    "scrollX": true
  });

});
