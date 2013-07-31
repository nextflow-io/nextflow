import nextflow.util.DxFile

params.genome= 'file-B7jPG5j0FqX1ByBjfj4005xV'
genome = new DxFile(id: params.genome, name:'genome')

params.annotation= "/usr/bin/annotation"
annotation = new File(params.annotation)

params.primary= "/usr/bin/primary"
primary = new File(params.primary)

task {
    echo true

    input genome
    output file1


    """
    cp ${genome} file1
    """

}

task {
  input file1
  output '*.txt'

  """
  echo 'This is the second task'
  cat ${file1} > prueba.txt
  """

}



