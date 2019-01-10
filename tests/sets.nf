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

x = Channel.from( ['a', 'file1'], ['b','file2'] )

process touch {

  input:
    set ( id, fileName ) from x
  output:
    set ( id, 'file*' ) into z


  /
  echo Creating $id
  touch $fileName
  /
}

process makeFiles {
  input:
    set( id, 'file_x' ) from z

  output:
    set( id, '*') into q mode flatten

  /
   cp file_x copy_$id
   touch beta_$id
  /

}

q.subscribe { println it }
