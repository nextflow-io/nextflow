process stubtest {
    debug true
    errorStrategy 'retry'

    output:
    path("*.txt")

    script:
    """
    echo "Not stubbing"
    touch script.txt
    """

    stub:
    if( task.attempt < 2 ) {
    """
    echo "Stubbing. Not creating file"
    """
    } else {
    """
    echo "Stubbing. Creating file"
    touch script.txt
    """
    }
}

workflow {
    main:
    stubtest()
}
