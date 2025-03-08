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
   channel.of('a'..'z') | foo | view
}