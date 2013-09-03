#!/usr/bin/env nextflow

params.in = "$HOME/sample.fa"

sequences = file(params.in)
SPLIT = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

task {
  input file:'input.fa', from: sequences
  output file:'seq_*', into: records

  """
  $SPLIT input.fa '%^>%' '/^>/' '{*}' -f seq_
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