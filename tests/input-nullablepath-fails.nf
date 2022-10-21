nextflow.enable.dsl=2
process test_process1 {
  input:
    val id
  output:
    path("output.txt", nullable:true)
  exec:
    println id
}

process test_process2 {
  input:
    path(file)
  output:
    val file
  exec:
    sleep 1000L
    println file
}

workflow {
    channel.of('foo') | test_process1 | test_process2 |  view()
}