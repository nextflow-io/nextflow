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
 * Collect all tuples with the same 'barcode'
 */

aggregation = algn_files.groupTuple()

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
