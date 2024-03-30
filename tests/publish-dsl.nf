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
nextflow.preview.topic = true

process align {
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
  output:
  path 'xxx'

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

  topic:
  foo.out >> 'foo'
}

output {
  directory 'results'

  'data' {
    from align.out[0], mode: 'copy'
    from align.out[1], mode: 'copy'
    from my_combine.out
    from 'foo', mode: 'link'
  }

  'data/more' {
    defaults mode: 'copy'
    from my_combine.out
  }
}
