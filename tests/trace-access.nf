process foo {
  memory { (trace != null) ? trace.memory * 2 : (1.GB) }
  errorStrategy 'retry'
  maxRetries 3
  input:
     val i
  output:
     stdout
  script:
  if( task.attempt <= 3 ){
  """
  echo mem: $task.memory
  exit 137
  """
  } else {
  """
  echo mem: $task.memory
  exit 0
  """
  }
}

workflow {
  foo(channel.of(1)).view() 
}
