workflow {
    foo | mix | collect | bar
}

process foo {
  publishDir 'foo'
  container 'debian:latest'
  output:
  file '*.fa'
  file 'hello.txt'
  file '*.{zip,html}'
  file '01_A(R{1,2}).fastq'
  file 'sample_(1 2).vcf'
  file '.alpha'

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
  debug true
  container 'debian:latest'
  input: 
  file '*'

  script:
  $/
  cat .alpha/hello.txt
  [ `cat * | grep -c ''` == 9 ] || false
  /$
}
