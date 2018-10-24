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
  * Run a process using a S3 file as input 
  */


s3file = file('s3://cbcrg-eu/ggal/ggal_1_48850000_49020000.bed.gff')
s3glob = Channel.fromFilePairs('s3://cbcrg-eu/ggal/*_{1,2}.fq')

process foo {
  echo true
  input:
  file(obj) from s3file

  """
  cat $obj | head
  """
}

process bar {
  tag "$pair"
  input:
  set pair, file(obj) from s3glob

  """
  cat $obj | head
  """

}