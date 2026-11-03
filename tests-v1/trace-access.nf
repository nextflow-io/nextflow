process foo {
  memory { task.attempt > 1 ? task.previousTrace.memory * 2 : (1.GB) }
  errorStrategy 'retry'
  maxRetries 3
  input:
     val i
  output:
     stdout
  script:
  if( task.attempt <= 3 ){
  """
  exit 137
  """
  } else {
  """
  echo 'mem: $task.memory (previous: $task.previousTrace.memory) (error: $task.previousException)'
  exit 0
  """
  }
}

workflow {
  foo(channel.of(1)).view() 
}
