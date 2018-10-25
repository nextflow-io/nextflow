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

def rule( file ) {
  if( file == 'file_1.txt' )
    return "alpha/$file"

  if( file == 'file_2.txt' )
    return null

  if( file == 'file_3.txt' )
     return "$PWD/results/gamma/$file"

}

process foo {
  publishDir path: 'results', saveAs: this.&rule

  input: each x from 1,2,3
  output: file '*.txt'
  """
  touch file_${x}.txt
  """

}
