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

x = 100
y = 200

process foo {
  input:
  file fastq from 'dummy'

  output:
  val 'Hello' into str_channel
  val "${fastq.baseName}-${x}.out" into exp_channel
  val x into x_channel
  val y into y_channel

  script:
  y = 'two hundred'
  """
  echo bar
  """
}

x_channel.println { "x: $it" }
y_channel.println { "y: $it" }
str_channel.println { "str: $it" }
exp_channel.println { "exp: $it" }
