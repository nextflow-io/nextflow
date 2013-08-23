#!/usr/bin/env nextflow

params.in = './sample.fa'

sequences = file(params.in)
SPLIT = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

task {
  input sequences
  output file:'seq_*', into: records
  
  """
  $SPLIT ${sequences} '%^>%' '/^>/' '{*}' -f seq_
  """

}

task {
    input records
    stdout reverse

    """
    cat $records | rev
    """
}

reverse.each { println it }