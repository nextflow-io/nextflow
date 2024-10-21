process foo {
  fair true
  input:
  val x 
  output: 
  tuple val(task.index), val(x)

  script:
  """
  sleep \$((RANDOM % 3))
  """
}


workflow {
  Channel.of('a'..'z') | foo | view
}