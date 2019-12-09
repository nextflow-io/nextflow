process foo {
    output:
    env FOO into ch 
    /FOO=Hello/
}

process bar {
    echo true
    input:
    env FOO from ch 
    'echo "bar says $FOO"'
}