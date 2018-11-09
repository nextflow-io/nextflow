#!/usr/bin/env nextflow
/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
 
_in1 = Channel.fromPath("$baseDir/data/file\\[a-b\\].txt")

process foo {

  input:
  file x from _in1

  output:
  file x into _out1
  file 'file-\\*.txt' into _out2
  file 'file-?.txt' glob false into _out3

  '''
  touch file-\\*.txt
  touch file-\\?.txt
  '''

}

_out1.println { "match: ${it.name}" }
_out2.println { "match: ${it.name}" }
_out3.println { "match: ${it.name}" }
