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
