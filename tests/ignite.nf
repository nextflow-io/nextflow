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

cheers = Channel.from 'Bonjour', 'Ciao', 'Hello', 'Hola', 'Γεια σου'

process sayHello {
  echo true

  input:
  val x from cheers

  """
  echo '$x world!'
  """
}

process printHello {
  input:
  val x from 1,2,3

  output:
  val z into results

  exec:
  println "Native: $x"
  z = x*x

}

results.println { "Square: $it" }
