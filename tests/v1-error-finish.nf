#!/usr/bin/env nextflow
// parser_v1: mixing script declarations and statements

process foo {
  debug true
  errorStrategy 'finish'
  input:  each x
  output: stdout

  script:
  if( x != 3 )
  """
    echo run_$x  
    sleep 5
  """
  else
  """
    exit 99
  """
}

process bar {
  input:  file 'x'

  script:
  '''
  cat x
  '''
}


workflow.onError {
  println "success: $workflow.success"
  println "exitStatus: $workflow.exitStatus"
}

workflow {
  foo([1,2,3]) | bar
}
