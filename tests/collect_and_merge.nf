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

/*
 * fake alignment step producing a BAM and BAI files
 */
process algn {
  echo true

  input:
  each barcode from 'alpha', 'gamma'
  each seq_id from 'one', 'two', 'three'

  output:
  set barcode, seq_id, file('bam'), file('bai') into algn_files

  """
  echo BAM $seq_id - $barcode > bam
  echo BAI $seq_id - $barcode > bai

  """

}

/*
 * aggregation is made by using a 'reduce' operator
 * followed by 'flatMpa'
 */

aggregation = algn_files
                .reduce([:]) { map, tuple ->    // 'map' is used to collect all values; 'tuple' is the record containing four items: barcode, seqid, bam file and bai file
                    def barcode = tuple[0]      // the first item is the 'barcode'
                    def group = map[barcode]    // get the aggregation for current 'barcode'
                    if( !group ) group = [ barcode, [], [], [] ]    // if new, create a new entry
                    group[1] << tuple[1]        // append 'seq_id' to the aggregation list
                    group[2] << tuple[2]        // append 'bam' file to the aggregation list
                    group[3] << tuple[3]        // append 'bai' file to the aggregation list
                    map[barcode] = group        // set back into the map
                    return map                  // return it so that it will be used in the next iteration
                }
                .flatMap { it.values() }        // tricky part: get the list of values of in the map, each value is the
                                                // aggregation build above
                                                // the 'flatMap' emits each of these aggregation list as a single item

                .map { it.collect {  it instanceof Collection ? it.sort() : it }   }

/*
 * Finally merge the BAMs and BAIs with the same 'barcode'
 */

process merge {

  echo true

  input:
  set barcode, seq_id, file(bam: 'bam?'), file(bai: 'bai?') from aggregation

  """
  echo barcode: $barcode
  echo seq_ids: $seq_id
  echo bam    : $bam
  echo bai    : $bai
  """

}
