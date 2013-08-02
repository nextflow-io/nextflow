import nextflow.util.DxFile

params.genome= 'file-B7jPG5j0FqX1ByBjfj4005xV'
genome = new DxFile(id: params.genome, name:'genome')
//genome= file(params.genome)

echo true

out = channel()
prueba = channel()


task {
  input genome
  output '*.txt': out
  output '*fa': prueba

  """
  echo 'This is the second task' >> out.txt
  cat $genome >> prueba.fa
  """

}


task {
  input out
  input prueba
  output reverse

  """
  cat $prueba | rev >> reverse
  cat $out >> reverse
  cat reverse
  """
}