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

seqs = channel.fromList(file("$baseDir/data/seqs/*.fastq"))

workflow {
  seqs | proc1
  seqs | proc2
  seqs | proc3
  seqs | proc4
  seqs | proc5
  seqs | proc6
  seqs | proc7
}

process proc1 {
  storeDir 'results'

  input:
  file fastq

  output:
  tuple val(baseName), file(outName) 

  script:
  baseName = fastq.baseName
  outName = baseName + '.proc1'
  """
  echo "blarg" > $outName
  """
}

process proc2 {
  storeDir 'results'

  input:
  file fastq

  output:
  tuple val(baseName), file(outName)

  script:
  baseName = fastq.baseName
  outName = baseName + '.proc2'
  """
  echo "blarg" > $outName
  """
}

process proc3 {
  storeDir 'results'

  input:
  file fastq

  output:
  tuple val(baseName), file(outName)

  script:
  baseName = fastq.baseName
  outName = baseName + '.proc3'
  """
  echo "blarg" > $outName
  """
}

process proc4 {
  storeDir 'results'

  input:
  file fastq

  output:
  tuple val(baseName), file(outName)

  script:
  baseName = fastq.baseName
  outName = baseName + '.proc4'
  """
  echo "blarg" > $outName
  """
}


process proc5 {
  storeDir 'results'

  input:
  file fastq

  output:
  tuple val(baseName), file(outName)

  script:
  baseName = fastq.baseName
  outName = baseName + '.proc5'
  """
  echo "blarg" > $outName
  """
}

process proc6 {
  storeDir 'results'

  input:
  file fastq

  output:
  tuple val(baseName), file(outName)

  script:
  baseName = fastq.baseName
  outName = baseName + '.proc6'
  """
  echo "blarg" > $outName
  """
}

process proc7 {
  storeDir 'results'

  input:
  file fastq

  output:
  tuple val(baseName), file(outName)

  script:
  baseName = fastq.baseName
  outName = baseName + '.proc7'
  """
  echo "blarg" > $outName
  """
}
