nextflow.preview.dsl=2 

process foo1 {
    echo true
    input: path x
    input: path y
    """
    echo "FOO1: ${x}; ${y}"
    """
}

process foo2 {
    echo true
    input: path x
    input: path y
    script:
    """
    echo "FOO2: ${x}; ${y}"
    """
}

process foo3 {
    echo true
    input: path x
    input: path y
    shell:
    '''
     echo "FOO3: !{x}; !{y}"
    '''
}

process foo4 {
    echo true
    input: path x
    input: path y
    script:
    template("$baseDir/task-escape-path-dsl2.sh")
}

workflow {
    f1 = file('file AA.txt')
    ch = Channel.fromPath(['file1.txt', 'file2.txt', 'fil BB.txt']).collect()
    foo1(f1,ch)
    foo2(f1,ch)
    foo3(f1,ch)
    foo4(f1,ch)
}