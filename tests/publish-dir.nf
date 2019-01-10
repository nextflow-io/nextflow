#!/usr/bin/env nextflow
/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
