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

def rule( file ) {
  if( file == 'file_1.txt' )
    return "alpha/$file"

  if( file == 'file_2.txt' )
    return null

  if( file == 'file_3.txt' )
     return "$PWD/results/gamma/$file"

}

process foo {
  publishDir path: 'results', saveAs: this.&rule

  input: each x from 1,2,3
  output: file '*.txt'
  """
  touch file_${x}.txt
  """

}
