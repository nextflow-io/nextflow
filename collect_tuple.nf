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
