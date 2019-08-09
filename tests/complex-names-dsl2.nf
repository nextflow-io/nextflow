#!/usr/bin/env nextflow
nextflow.preview.dsl=2

process foo {
  publishDir 'foo', mode: 'copy'
  container 'debian:latest'
  output:
  path '*.fa'
  path 'hello.txt'
  path '*.{zip,html}'
  path '01_A(R{1,2}).fastq'
  path 'sample_(1 2).vcf'
  path '.alpha'

  script:
  $/
  echo A > hello.txt
  echo B > sample.zip 
  echo C > sample.html
  echo D > 01_A\(R1\).fastq
  echo E > 01_A\(R2\).fastq
  echo F > sample_\(1\ 2\).vcf
  echo 1 > f1.fa
  echo 2 > f2.fa
  echo 3 > f3.fa
  mkdir .alpha
  echo "Hello world!" > .alpha/hello.txt
  /$
}

process bar {
  echo true
  container 'debian:latest'
  input: 
  path '*'

  script:
  $/
  cat .alpha/hello.txt
  [ `cat * | grep -c ''` == 9 ] || false
  /$
}

/*
 * main flow
 */

workflow {
    foo | mix | collect | bar
}
