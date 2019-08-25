#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
 * fake alignment step producing a BAM and BAI files
 */
process algn {
  echo true

  input:
  each barcode
  each seq_id

  output:
  set barcode, seq_id, path('bam'), path('bai')

  """
  echo BAM $seq_id - $barcode > bam
  echo BAI $seq_id - $barcode > bai
  """
}


/*
 * Finally merge the BAMs and BAIs with the same 'barcode'
 */

process merge {
  echo true

  input:
  set barcode, seq_id, path(bam, stageAs:'bam?'), path(bai, stageAs:'bai?')

  """
  echo barcode: $barcode
  echo seq_ids: $seq_id
  echo bam    : $bam
  echo bai    : $bai
  """
}

/*
 * main flow
 */

workflow {
    algn( ['alpha', 'gamma'], ['one', 'two', 'three'] ) \
      | groupTuple \
      | merge

}
