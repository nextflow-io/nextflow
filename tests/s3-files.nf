#!/usr/bin/env nextflow
/*
 * Copyright 2020, Seqera Labs
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
  * Run a process using a S3 file as input 
  */


s3file = file('s3://nextflow-ci/nf-test-data/transcriptome.fa')
s3glob = Channel.fromFilePairs('s3://nextflow-ci/nf-test-data/*_{1,2}.fq')

process foo {
  echo true
  input:
  file(obj) from s3file

  """
  cat $obj | head
  """
}

process bar {
  tag "$pair"
  input:
  set pair, file(obj) from s3glob

  """
  cat $obj | head
  """

}
