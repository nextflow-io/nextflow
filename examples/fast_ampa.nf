#!/bin/env nextflow

if( !args ) {
  exit "Please specify the input file as a cmd line argument"
}

fileName = args[0]
fastaFile = new File(fileName)

seq = channel()
fastaFile.chunkFasta { seq << it }

ampaOut = channel()
task ('ampa') {

    //  defines the Input and Output
    //echo true
    input '-':seq
    output '-':ampaOut
    threads 10

    // The BASH script to be executed - for each - sequence
    """
    set -e
    cat - > input.file && AMPA-BIGTABLE.pl -in=input.file -noplot -rf=result -df=data
    cat input.file | grep '>' >> /dev/stdout
    cat result | grep '#' >> /dev/stdout
    """

}

def result = new File('bigampa.txt')
println "Saving result at: ${result.absoluteFile}"

ampaOut.eachWithIndex { str, index ->

  def lines = str.trim().split('\n')
  if( lines.size() != 2 ) {
      println "ERROR > Invalid AMPA ($index) result:\n${str}\n\n"
      return
  }

  def id = getIDs(lines[0])
  def val = getValues(lines[1])

  result << "${id[0]}\t${id[1]}\t${val[0]}\t${val[1]}\t${val[2]}\t${val[3]}\n"
}



def getIDs( line ) {

  def matcher = line =~ />(\S+).+gene:(\S+).*/
  if( matcher.matches() ) {
    def seqId = matcher[0][1]
    def geneId = matcher[0][2]
    return [seqId, geneId]
  }
  else {
    println "Bad ID line: $line"
    return []
  }

}


/*
 *  return the values in the following order 
 *  - stretches 
 *  - protLen: 
 *  - ampLen: 
 *  - propensity
 */
def getValues(result) {

 def rm = result =~ /# This protein has (\d+) bactericidal stretches and it has (\d+) amino acids. AMP length: (\d+) Best AMP Propensity: ([0-9\.]+)/

  if( rm.matches() ) {
    return [rm[0][1], rm[0][2], rm[0][3], rm[0][4]]
  }
  else {
    println "Bad result line: $result"
    return []
  }
}


