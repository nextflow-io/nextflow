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


process hola {
    debug true
    input:
    val x
    each y
    each z

    """
    echo 'x: $x; y: $y; z: $z'
    """

}

process foo {
    debug true

    input:
    each v

    """
    echo foo $v
    """
}

workflow {
  def list1 = channel.of(1,2)
  def list2 = channel.of('Hola', 'Ciao')
  def list3 = channel.of('alpha','beta','delta')
  def list4 = channel.of(["a","b"],["c","d"])

  hola(list1, list2, list3)

  foo(list4)
}
