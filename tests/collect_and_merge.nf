#!/usr/bin/env nextflow
/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
