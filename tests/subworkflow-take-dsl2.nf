nextflow.preview.dsl=2

process foo {
    input: val x 
    output: val x
    /echo true/
}

process bar { 
    input: val x 
    output: val x
    /echo true/
}

workflow flow1 {
    take: data
    main:
        foo(data)
        bar(foo.out)
    emit:
        bar.out
}

workflow flow2 {
    take: data
    main:
        foo(data)
        bar(foo.out)
    emit:
        bar.out
}

workflow {
    flow1('foo')
    flow2(flow1.out)
}