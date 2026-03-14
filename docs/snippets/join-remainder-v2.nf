nextflow.preview.types = true

workflow {
    left  = channel.of( record(id: 'X', a: 1), record(id: 'Y', a: 2), record(id: 'Z', a: 3), record(id: 'P', a: 7) )
    right = channel.of( record(id: 'Z', b: 6), record(id: 'Y', b: 5), record(id: 'X', b: 4), record(id: 'Q', b: 8) )

    left.join(right, by: 'id', remainder: true).view()
}