nextflow.enable.dsl=1

process foo {
    output:
    env FOO into ch 
    /FOO=Hello/
}

process bar {
    debug true
    input:
    env FOO from ch 
    'echo "bar says $FOO"'
}
