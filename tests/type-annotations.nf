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

workflow {
  words_ch = channel.of('one', 'two', 'three', 'four')
  counts_vl = COUNT(words_ch)
  counts_vl.view { counts ->
    def even = counts.findAll { n -> isEven(n) }.size()
    println "counts: $counts ($even are even)"
  }
}

workflow COUNT {
  take:
  words: Channel<String>

  main:
  counts = words.map { word -> word.length() }.collect()

  emit:
  counts: Value<Integer> = counts
}

def isEven(n: Integer) -> Boolean {
  return n % 2 == 0
}
