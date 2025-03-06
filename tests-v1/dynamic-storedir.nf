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

params.prefix = 'my'

process foo {
  storeDir "cache/$x"

  input:
  each x
  file 'result.txt'

  output:
  tuple val(x), file('result.txt')

  """
  echo World >> result.txt
  """

}

workflow {
  def data = 'Hello\n'
  def list = ['alpha', 'delta', 'gamma', 'omega']
  foo(list, data) | subscribe { code, file ->
                      println "~ Result ${file}"
                      file.copyTo("my_${code}.txt")
                    }

}
