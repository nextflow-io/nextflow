#!/usr/bin/env nextflow
/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

input = Channel.from('alpha','beta','delta')

process align {
    publishDir 'data', mode: 'copy'

    input:
    val(x) from input

    output:
    file("*.bam") into bam
    file("${x}.bai") into bai

    """
    echo ${x} > ${x}.bam
    echo ${x} | rev > ${x}.bai
    """
}

process combine {
    publishDir 'data'
    publishDir 'more/data', mode: 'copy'

    input:
    file(bamfile) from bam.toSortedList { it.name }
    file(baifile) from bai.toSortedList { it.name }

    output:
    file 'result.txt' into result

    """
    cat $bamfile > result.txt
    cat $baifile >> result.txt
    """
}

result.subscribe { println it.text }

process foo {
  publishDir 'data', mode: 'link'
  output:
  file 'xxx'

  '''
  mkdir xxx
  touch xxx/A
  touch xxx/B
  touch xxx/C
  '''
}
