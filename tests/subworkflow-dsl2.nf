#!/bin/bash nextflow
nextflow.preview.dsl=2

process foo {
    output: stdout()
    shell: "echo Hello"
}

process bar {
    input: file "foo"
    output: stdout()
    shell:
    "rev foo"
}

process baz {
    input: file "foo"
    output: stdout()
    shell: "tr '[:lower:]' '[:upper:]' < foo"
}

workflow flow1 {
    main:
    foo | bar | collectFile | set { result }
    emit: 
    result
}

workflow flow2 {
    main:
    foo | baz | collectFile | set { result }
    emit: 
    result
}

workflow test1 {
    flow1()
    flow2()
    ch1 = flow1.out.result
    ch2 = flow2.out.result
    emit: ch1.mix(ch2).collectFile(name:"$PWD/test1.txt")
}

workflow test2 {
    emit: ( flow1 & flow2 ) | mix | collectFile(name:"$PWD/test2.txt")
}

workflow {
    test1().view()
    test2().view()
}
