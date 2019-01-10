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
 
process foo {
  output:
  file 'x' into A,B,C
  val 2 into baz

  """
  echo Ciao > x
  """
}

baz.println { "Hello $it" }

process bar1 {
  echo true
  input:
  file x from A

  """
  printf "\$(cat $x) A"
  """
}

process bar2 {
  echo true
  input:
  file x from B

  """
  printf "\$(cat $x) B"
  """
}

process bar3 {
  echo true
  input:
  file x from C

  """
  printf "\$(cat $x) C"
  """
}
