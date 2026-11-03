#!/usr/bin/env nextflow
/*
 * Copyright 2013-2024, Seqera Labs
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


process align {
    publishDir 'data', mode: 'copy'

    input:
    val(x)

    output:
    path("*.bam")
    path("${x}.bai")

    """
    echo ${x} > ${x}.bam
    echo ${x} | rev > ${x}.bai
    """
}

process my_combine {
    publishDir 'data'
    publishDir 'more/data', mode: 'copy'

    input:
    path(bamfile)
    path(baifile)

    output:
    path 'result.txt'

    """
    cat $bamfile > result.txt
    cat $baifile >> result.txt
    """
}

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

workflow {
  def input = Channel.of('alpha','beta','delta')
  align(input)

  def bam = align.out[0].toSortedList { it.name }
  def bai = align.out[1].toSortedList { it.name }
  my_combine( bam, bai )
  my_combine.out.view{ it.text }

  foo()
}
