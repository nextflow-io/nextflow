#!/usr/bin/env nextflow

process foo {
  output:
    path('one.txt', arity: '1')
    path('pair_*.txt', arity: '2')
    path('many_*.txt', arity: '1..*')
  script:
    """
    echo 'one' > one.txt
    echo 'pair_1' > pair_1.txt
    echo 'pair_2' > pair_2.txt
    echo 'many_1' > many_1.txt
    echo 'many_2' > many_2.txt
    echo 'many_3' > many_3.txt
    """
}

process bar {
  input:
    path('one.txt', arity: '1')
    path('pair_*.txt', arity: '2')
    path('many_*.txt', arity: '1..*')
  script:
    """
    cat one.txt
    cat pair_*.txt
    cat many_*.txt
    """
}

workflow {
    foo | bar
}
