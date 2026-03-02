#!/usr/bin/env nextflow

nextflow.preview.types = true

/*
 * fake alignment step producing a BAM and BAI files
 */
process align {
  debug true

  input:
  in: Record {
    barcode: String
    seq_id: String
  }

  output:
  record(
    barcode: barcode,
    seq_id: seq_id,
    bam: file('bam'),
    bai: file('bai')
  )

  script:
  """
  echo BAM ${in.seq_id} - ${in.barcode} > bam
  echo BAI ${in.seq_id} - ${in.barcode} > bai
  """
}


/*
 * Finally merge the BAMs and BAIs with the same 'barcode'
 */
process merge {
  debug true

  input:
  group: Record {
    barcode: String
    seq_ids: Bag<String>
    bam: Bag<Path>
    bai: Bag<Path>
  }

  stage:
  stageAs 'bam?', group.bam
  stageAs 'bai?', group.bai

  script:
  """
  echo barcode: ${group.barcode}
  echo seq_ids: ${group.seq_ids.join(' ')}
  echo bam    : ${group.bam.join(' ')}
  echo bai    : ${group.bai.join(' ')}
  """
}


/*
 * main flow
 */
workflow {
  ch_barcode = channel.of('alpha', 'gamma')
  ch_seq = channel.of('one', 'two', 'three')
  ch_inputs = ch_barcode
    .cross(ch_seq)
    .map { barcode, seq_id -> record(barcode: barcode, seq_id: seq_id) }

  ch_aligned = align( ch_inputs )
  ch_grouped = ch_aligned.groupTuple { r -> tuple(r.barcode, r) }
  merge( ch_grouped )
}
