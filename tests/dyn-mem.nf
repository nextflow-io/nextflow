
process foo {
  memory { x.size() < 10.B  ? 100.MB : 200.MB }
  
  input: 
  file x from Channel.fromPath(['.small.txt','.big.txt'])
  
  script:
  """
  echo task=$x mem=$task.memory 
  """
}