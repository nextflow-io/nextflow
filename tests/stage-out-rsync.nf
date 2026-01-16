
process test {
    tag "$number"
    storeDir 'cache/test'
    stageOutMode 'rsync'

    input:
    val number

    output:
    path('*.txt')

    script:
    """
    echo "This is an example process with number: ${number}" > ${number}.txt
    """
}

workflow {
    numbers = channel.of(1, 2, 3)
    test(numbers).view()
}
