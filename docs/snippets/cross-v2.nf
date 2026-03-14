nextflow.preview.types = true

workflow {
    numbers = channel.of(1, 2, 3)
    words = channel.of('hello', 'ciao')

    numbers.cross(words).view()
}