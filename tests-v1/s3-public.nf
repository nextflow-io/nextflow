process foo {
  input:
  path x
  script:
  """
  cat $x
  """
}

workflow {
  def f0 = file('s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/README.txt')
  foo(f0)
}