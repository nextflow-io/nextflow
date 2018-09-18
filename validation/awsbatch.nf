process foo {
  publishDir 'foo'
  container 'debian:latest'
  output:
  file '*.fa' into ch1
  file 'hello.txt'  into ch2 
  file '*.{zip,html}' into ch3  
  file '01_A(R{1,2}).fastq' into ch4 
  file 'sample_(1 2).vcf' into ch5
  file '.alpha' into ch6

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
  file '*' from ch1.mix(ch2,ch3,ch4,ch5,ch6).collect()

  script:
  $/
  cat .alpha/hello.txt
  [ `cat * | grep -c ''` == 9 ] || false
  /$
}
