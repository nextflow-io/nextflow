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

process touch {
  input:
    tuple val(id), val(fileName)
  output:
    tuple val(id), path('file*')

  /
  echo Creating $id
  touch $fileName
  /
}

process makeFiles {
  input:
    tuple val(id), path('file_x')

  output:
    tuple val(id), path('*')

  /
   cp file_x copy_$id
   touch beta_$id
  /
}


workflow {
  def x = Channel.from( ['a', 'file1'], ['b','file2'] )
  touch(x)
  makeFiles(touch.out)
  makeFiles.out.view()
}
