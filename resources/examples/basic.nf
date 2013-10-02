#!/usr/bin/env nextflow

params.in = "$HOME/sample.fa"

sequences = file(params.in)
records = channel()
SPLIT = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

task {
  input sequences
  output 'seq_*': records
  
  """
  $SPLIT $sequences '%^>%' '/^>/' '{*}' -f seq_
  """

}


reverse = channel()

task {
    input records
    output '-': reverse

    """
    cat $records | rev
    """
}


reverse.each { println it }