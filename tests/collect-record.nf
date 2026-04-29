#!/usr/bin/env nextflow

nextflow.enable.types = true

/*
 * fake alignment step producing a BAM and BAI files
 */
process align {
    debug true

    input:
    record(
        barcode: String,
        seq_id: String
    )

    output:
    record(
        barcode: barcode,
        seq_id: seq_id,
        bam: file('bam'),
        bai: file('bai')
    )

    script:
    """
    echo BAM ${seq_id} - ${barcode} > bam
    echo BAI ${seq_id} - ${barcode} > bai
    """
}

/*
 * Finally merge the BAMs and BAIs with the same 'barcode'
 */
process merge {
    debug true

    input:
    tuple(barcode: String, samples: Bag<Record>)

    stage:
    stageAs samples*.bam, 'bam?'
    stageAs samples*.bai, 'bai?'

    script:
    """
    echo barcode: ${barcode}
    echo seq_ids: ${samples*.seq_id.join(' ')}
    echo bam    : ${samples*.bam.join(' ')}
    echo bai    : ${samples*.bai.join(' ')}
    """
}

/*
 * main flow
 */
workflow {
    ch_barcode = channel.of('alpha', 'gamma')
    ch_seq = channel.of('one', 'two', 'three')
    ch_inputs = ch_barcode
        .combine(ch_seq)
        .map { barcode, seq_id -> record(barcode: barcode, seq_id: seq_id) }

    ch_aligned = align( ch_inputs )
        .map { r -> tuple(r.barcode, r) }
        .groupBy()

    merge( ch_aligned )
}
