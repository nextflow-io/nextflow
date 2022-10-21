nextflow.enable.dsl=2

process test_process1 {
  input:
    val id
  output:
    path("output.txt")
  exec:
    println 'hi'
}
workflow {
    test_process1('foo').out
}