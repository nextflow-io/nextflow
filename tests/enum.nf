
enum MyEnum {
    A,
    B
}

workflow {
    ch = channel.of( tuple([:], MyEnum.B) )
    SomeTask(ch)
}

process SomeTask {
    input:
    tuple val(meta), val(myEnum)

    output:
    tuple val(meta), path('out.txt')

    script:
    """
    echo "hello $myEnum" > out.txt
    """
}
