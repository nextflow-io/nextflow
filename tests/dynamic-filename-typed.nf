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

nextflow.preview.types = true

params {
  input: Path
  prefix: String = 'my'
  names: List<String> = ['alpha', 'delta', 'gamma', 'omega']
}

process foo {
  stageInMode 'copy'

  input:
  (name, txt): Tuple<String, Path>
  prefix: String

  stage:
  stageAs "${prefix}_${name}.txt", txt

  output:
  file("${prefix}_${name}.txt")

  script:
  """
  echo World >> ${prefix}_${name}.txt
  """
}

workflow {
  ch_names = channel.fromList(params.names)
  ch_foo = foo( ch_names.cross(params.input) )
  ch_foo.subscribe { it ->
    println "~ Saving ${it.name}"
    it.copyTo('.')
  }
}
