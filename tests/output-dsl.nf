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
nextflow.preview.output = true

params.save_foo = true

process align {
  input:
  val(x)

  output:
  path("*.bam")
  path("${x}.bai")

  script:
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

  script:
  """
  cat $bamfile > result.txt
  cat $baifile >> result.txt
  """
}

process foo {
  output:
  path 'xxx'

  script:
  '''
  mkdir xxx
  touch xxx/A
  touch xxx/B
  touch xxx/C
  '''
}

workflow {
  main:
  input = Channel.of('alpha','beta','delta')
  align(input)

  bams = align.out[0].toSortedList { bam -> bam.name }
  bais = align.out[1].toSortedList { bai -> bai.name }
  my_combine( bams, bais )
  my_combine.out.view { it -> it.text }

  foo()

  publish:
  align.out       >> 'data'
  my_combine.out  >> 'more/data'
  foo.out         >> (params.save_foo ? 'data' : null)
}

output {
  'data' {
    path { val -> { file -> file } }
    index {
      path 'index.csv'
      mapper { val -> [filename: val] }
      header true
      sep ','
    }
  }
}
