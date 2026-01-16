#!/usr/bin/env nextflow

nextflow.preview.types = true

/*
 * fake alignment step producing a BAM and BAI files
 */
process align {
  debug true

  input:
  (barcode, seq_id): Tuple<String, String>

  output:
  tuple(barcode, seq_id, file('bam'), file('bai'))

  script:
  """
  echo BAM $seq_id - $barcode > bam
  echo BAI $seq_id - $barcode > bai
  """
}


/*
 * Finally merge the BAMs and BAIs with the same 'barcode'
 */
process merge {
  debug true

  input:
  (barcode, seq_ids, bam, bai): Tuple<String, Bag<String>, Bag<Path>, Bag<Path>>

  stage:
  stageAs 'bam?', bam
  stageAs 'bai?', bai

  script:
  """
  echo barcode: $barcode
  echo seq_ids: ${seq_ids.join(' ')}
  echo bam    : ${bam.join(' ')}
  echo bai    : ${bai.join(' ')}
  """
}


/*
 * main flow
 */
workflow {
  ch_barcode = channel.of('alpha', 'gamma')
  ch_seq = channel.of('one', 'two', 'three')

  align( ch_barcode.combine(ch_seq) )
    | groupTuple
    | merge
}
