// JavaScript used to power the Nextflow Report Template output.
window.data_byprocess = {};

/* helper functions that takes an array of numbers each of each
   is a integer representing a number of bytes and normalise to base 2 scale */
function norm_mem( list ) {
  if( list == null ) return null;
  var result = new Array(list.length);
  for( i=0; i<list.length; i++ ) {
    var value = list[i];
    var x = Math.floor(Math.log10(value) / Math.log10(1024));
    if( x == 0 )
      value = value/1.024;
    else {
      for( j=0; j<x; j++ )
        value = value / 1.024;
    }
    result[i] = Math.round(value);
  }
  return result;
}

$(function() {
  // Script block clicked
  $('#tasks_table').on('click', '.script_block', function(e){
    e.preventDefault();
    $(this).toggleClass('short');
  });

  $(function() {
    $('[data-toggle="tooltip"]').tooltip()
  })

  // Completed date from now
  var completed_date = moment( $('#workflow_complete').text(), "ddd MMM DD HH:mm:ss .* YYYY" );
  if(completed_date.isValid()){
    $('#completed_fromnow').html('completed ' + completed_date.fromNow() + ', ');
  }

  // Collect metrics by process
  for(let i in window.data.summary){
    let metrics = window.data.summary[i];
    let proc = metrics.process;
    
    if(!window.data_byprocess.hasOwnProperty(proc)){
      window.data_byprocess[proc] = {};
    }

    for (let key in metrics) {
      if (metrics[key] != null) {
        window.data_byprocess[proc][key] = [];
        if( metrics[key].min == metrics[key].max ) {
            // min equals max ==> show just a value
            window.data_byprocess[proc][key].push(metrics[key].min);
        }
        else {
            // otherwise show all values
            window.data_byprocess[proc][key].push(metrics[key].min);
            window.data_byprocess[proc][key].push(metrics[key].q1);
            window.data_byprocess[proc][key].push(metrics[key].q1);
            window.data_byprocess[proc][key].push(metrics[key].q2);
            window.data_byprocess[proc][key].push(metrics[key].q3);
            window.data_byprocess[proc][key].push(metrics[key].q3);
            window.data_byprocess[proc][key].push(metrics[key].max);
        }
        if (key == "time") {
            window.data_byprocess[proc][key] = window.data_byprocess[proc][key].map(function(d,i){
            return moment.duration(d).asMinutes().toFixed(1);
          });
        }
      }
    }
  }

  // Plot histograms of resource usage
  var cpu_raw_data = [];
  var cpu_usage_data = [];
  var mem_raw_data = [];
  var mem_usage_data = [];
  var vmem_raw_data = [];
  var time_raw_data = [];
  var time_usage_data = [];
  var reads_raw_data = [];
  var writes_raw_data = [];
  for(var pname in window.data_byprocess){
    if( !window.data_byprocess.hasOwnProperty(pname) )
        continue;
    var smry = window.data_byprocess[pname];
    cpu_raw_data.push({y: smry.cpu, name: pname, type:'box', boxmean: true, boxpoints: false});
    cpu_usage_data.push({y: smry.cpuUsage, name: pname, type:'box', boxmean: true, boxpoints: false});
    mem_raw_data.push({y: norm_mem(smry.mem), name: pname, type:'box', boxmean: true, boxpoints: false});
    mem_usage_data.push({y: smry.memUsage, name: pname, type:'box', boxmean: true, boxpoints: false});
    vmem_raw_data.push({y: norm_mem(smry.vmem), name: pname, type:'box', boxmean: true, boxpoints: false});
    time_raw_data.push({y: smry.time, name: pname, type:'box', boxmean: true, boxpoints: false});
    time_usage_data.push({y: smry.timeUsage, name: pname, type:'box', boxmean: true, boxpoints: false});
    reads_raw_data.push({y: norm_mem(smry.reads), name: pname, type:'box', boxmean: true, boxpoints: false});
    writes_raw_data.push({y: norm_mem(smry.writes), name: pname, type:'box', boxmean: true, boxpoints: false});
  }

  Plotly.newPlot('cpuplot', cpu_raw_data, { title: 'CPU Usage', yaxis: {title: '% single core CPU usage', tickformat: '.1f', rangemode: 'tozero'} });
  Plotly.newPlot('memplot', mem_raw_data, { title: 'Physical Memory Usage', yaxis: {title: 'Memory', tickformat: '.4s', rangemode: 'tozero'} });
  Plotly.newPlot('timeplot', time_raw_data, { title: 'Task execution real-time', yaxis: {title: 'Execution time (minutes)', tickformat: '.1f', rangemode: 'tozero'} });
  Plotly.newPlot('readplot', reads_raw_data, { title: 'Number of bytes read', yaxis: {title: 'Read bytes', tickformat: '.4s', rangemode: 'tozero'} });

  // Only plot tabbed plots when shown
  $('#pctcpuplot_tablink').on('shown.bs.tab', function (e) {
    if($('#pctcpuplot').is(':empty')){
      Plotly.newPlot('pctcpuplot', cpu_usage_data, { title: '% Requested CPU Used', yaxis: {title: '% Allocated CPUs Used', tickformat: '.1f', rangemode: 'tozero'} });
    }
  });
  $('#pctmemplot_tablink').on('shown.bs.tab', function (e) {
    if($('#pctmemplot').is(':empty')){
        Plotly.newPlot('pctmemplot', mem_usage_data, { title: '% Requested Physical Memory Used', yaxis: {title: '% Memory', tickformat: '.1f', rangemode: 'tozero'} });
    }
  });
  $('#vmemplot_tablink').on('shown.bs.tab', function (e) {
    if($('#vmemplot').is(':empty')){
        Plotly.newPlot('vmemplot', vmem_raw_data, { title: 'Virtual Memory Usage', yaxis: {title: 'Memory', tickformat: '.4s', rangemode: 'tozero'} });
    }
  });
  $('#pcttimeplot_tablink').on('shown.bs.tab', function (e) {
    if($('#pcttimeplot').is(':empty')){
        Plotly.newPlot('pcttimeplot', time_usage_data, { title: '% Requested Time Used', yaxis: {title: '% Allocated Time Used', tickformat: '.1f', rangemode: 'tozero'} });
    }
  });
  $('#writeplot_tablink').on('shown.bs.tab', function (e) {
    if($('#writeplot').is(':empty')){
        Plotly.newPlot('writeplot', writes_raw_data, { title: 'Number of bytes written', yaxis: {title: 'Written bytes', tickformat: '.4s', rangemode: 'tozero'}});
    }
  });

  // Humanize duration
  function humanize(duration){
    if (duration.days() > 0) {
      return duration.days() + "d " + duration.hours() + "h"
    }
    if (duration.hours() > 0) {
      return duration.hours() + "h " + duration.minutes() + "m"
    }
    if (duration.minutes() > 0) {
      return duration.minutes() + "m " + duration.seconds() + "s"
    }
    return duration.asSeconds().toFixed(1) + "s"
  }

  // Build the trace table
  function make_duration(ms){
    if($('#nf-table-humanreadable').val() == 'false'){
      return ms;
    }
    if (ms == '-' || ms == 0){
      return ms;
    }
    return humanize(moment.duration( parseInt(ms) ));
  }
  function make_date(ms){
    if($('#nf-table-humanreadable').val() == 'false'){
      return ms;
    }
    if (ms == '-' || ms == 0){
      return ms;
    }
    return moment( parseInt(ms) ).format();
  }
  function make_memory(bytes){
    if($('#nf-table-humanreadable').val() == 'false'){
      return bytes;
    }
    if (bytes == '-' || bytes == 0){
      return bytes;
    }
    // https://stackoverflow.com/a/14919494
    var thresh = 1024;
    if(Math.abs(bytes) < thresh) {
      return bytes + ' B';
    }
    var units = ['kB','MB','GB','TB','PB','EB','ZB','YB'];
    var u = -1;
    do {
        bytes /= thresh;
        ++u;
    } while(Math.abs(bytes) >= thresh && u < units.length - 1);
    return bytes.toFixed(3)+' '+units[u];
  }
  function make_tasks_table(){
    // reset
      if ( $.fn.dataTable.isDataTable( '#tasks_table' ) ) {
        $('#tasks_table').DataTable().destroy();
      }
      var table = $('#tasks_table').DataTable({
        data: window.data.trace,
        columns: [
          { title: 'task_id', data: 'task_id' },
          { title: 'process', data: 'process' },
          { title: 'tag', data: 'tag' },
          { title: 'status', data: 'status', render: function(data, type, row){
              var s = {
                COMPLETED: 'success',
                CACHED: 'secondary',
                ABORTED: 'danger',
                FAILED: 'danger'
              }
              return '<span class="badge badge-'+s[data]+'">'+data+'</span>';
            }
          },
          { title: 'hash', data: 'hash', render:  function(data, type, row){
              var script = '';
              var lines = data.split("\n");
              var ws_re = /^(\s+)/g;
              var flws_match = ws_re.exec(lines[1]);
              if(flws_match == null){
                script = data;
              } else {
                for(var j=0; j<lines.length; j++){
                  script += lines[j].replace(new RegExp('^'+flws_match[1]), '').replace(/\s+$/,'') + "\n";
                }
              }
              return '<code>'+script+'</code>';
            }
          },
          { title: 'allocated cpus', data: 'cpus' },
          { title: '%cpu', data: '%cpu' },
          { title: 'allocated memory', data: 'memory', render: make_memory },
          { title: '%mem', data: '%mem' },
          { title: 'vmem', data: 'vmem', render: make_memory },
          { title: 'rss', data: 'rss', render: make_memory },
          { title: 'peak_vmem', data: 'peak_vmem', render: make_memory },
          { title: 'peak_rss', data: 'peak_rss', render: make_memory },
          { title: 'allocated time', data: 'time', render: make_duration },
          { title: 'duration', data: 'duration', render: make_duration },
          { title: 'realtime', data: 'realtime', render: make_duration },
          { title: 'script', data: 'script', render: function(data) {
              return '<pre class="script_block short"><code>' + data.trim() + '</code></pre>';
            }
          },
          { title: 'exit', data: 'exit' },
          { title: 'submit', data: 'submit', render: make_date },
          { title: 'start', data: 'start', render: make_date },
          { title: 'complete', data: 'complete', render: make_date },
          { title: 'rchar', data: 'rchar', render: make_memory },
          { title: 'wchar', data: 'wchar', render: make_memory },
          { title: 'syscr', data: 'syscr', render: make_memory },
          { title: 'syscw', data: 'syscw', render: make_memory },
          { title: 'read_bytes', data: 'read_bytes', render: make_memory },
          { title: 'write_bytes', data: 'write_bytes', render: make_memory },
          { title: 'native_id', data: 'native_id' },
          { title: 'name', data: 'name' },
          { title: 'module', data: 'module' },
          { title: 'container', data: 'container', render: function(data) {
              return '<samp>'+data+'</samp>';
            }
          },
          { title: 'disk', data: 'disk' },
          { title: 'attempt', data: 'attempt' },
          { title: 'scratch', data: 'scratch', render: function(data) {
              return '<samp>'+data+'</samp>';
            }
          },
          { title: 'workdir', data: 'workdir', render: function(data) {
              return '<samp>'+data+'</samp>';
            }
          }
        ],
        "deferRender": true,
        "lengthMenu": [[25, 50, 100, -1], [25, 50, 100, "All"]],
        "scrollX": true,
        "colReorder": true,
        "columnDefs": [
          { className: "id", "targets": [ 0,1,2,3 ] },
          { className: "meta", "targets": [ 4,13,16,17,18,19,20,27,28,29,30,31,32,33,34 ] },
          { className: "metrics", "targets": [ 5,6,7,8,9,10,11,12,14,15,21,22,23,24,25,26 ] }
        ],
        "buttons": [
          {
            extend: 'colvisGroup',
            text: 'Metrics',
            show: [ '.id', '.metrics' ],
            hide: [ '.meta' ],
          },
          {
            extend: 'colvisGroup',
            text: 'Metadata',
            show: [ '.id', '.meta'],
            hide: [ '.metrics' ],
          },
          {
            extend: 'colvisGroup',
            text: 'All',
            show: ':hidden',
          },
        ]
      });

      // Insert column filter button group
      table.buttons().container()
         .prependTo( $('#tasks_table_filter') );

      // Column filter button group onClick event to highlight active filter
      $('.buttons-colvisGroup').click(function(){
        var def = 'btn-secondary';
        var sel = 'btn-primary';
        $('.buttons-colvisGroup').removeClass(sel).addClass(def);
        $(this).addClass(sel).removeClass(def);
      });

      // Default filter highlight
      $(".buttons-colvisGroup:contains('All')").click();
    }

  if( window.data.trace==null ) {
      // nascondere
      $('#table-container').remove()
  }
  else {
      $('#no-table-container').remove()
      // Dropdown changed about raw / human readable values in table
      $('#nf-table-humanreadable').change(function(){
        make_tasks_table();
      });
      // Make the table on page load
      make_tasks_table();
  }


});
