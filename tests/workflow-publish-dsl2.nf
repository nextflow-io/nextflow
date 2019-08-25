nextflow.preview.dsl=2 

process foo {
  output: 
  file 'hello.txt' 
  /echo Hello world > hello.txt/
}

process bar {
  output: 
  file 'one.txt' 
  file 'two.txt'
  /
    echo Hello > one.txt
    echo world > two.txt
  /
}


workflow {
  main: 
    foo()
    bar()
  emit: 
    x = foo.out

  publish: 
    x to: 'results'
    bar.out to: 'results'
}
