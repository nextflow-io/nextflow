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

/*
 * @author Simon Ye
 */

seqs = file("$baseDir/data/seqs/*.fastq")

process proc1 {
  storeDir 'results'

  input:
  file fastq from seqs

  output:
  set val(baseName), file(outName) into proc1

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
  file fastq from seqs

  output:
  set val(baseName), file(outName) into proc2

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
  file fastq from seqs

  output:
  set val(baseName), file(outName) into proc3

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
  file fastq from seqs

  output:
  set val(baseName), file(outName) into proc4

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
  file fastq from seqs

  output:
  set val(baseName), file(outName) into proc5

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
  file fastq from seqs

  output:
  set val(baseName), file(outName) into proc6

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
  file fastq from seqs

  output:
  set val(baseName), file(outName) into proc7

  script:
  baseName = fastq.baseName
  outName = baseName + '.proc7'
  """
  echo "blarg" > $outName
  """
}
