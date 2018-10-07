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
 
_in1 = Channel.fromPath("$baseDir/data/file\\[a-b\\].txt")

process foo {

  input:
  file x from _in1

  output:
  file x into _out1
  file 'file-\\*.txt' into _out2
  file 'file-?.txt' glob false into _out3

  '''
  touch file-\\*.txt
  touch file-\\?.txt
  '''

}

_out1.println { "match: ${it.name}" }
_out2.println { "match: ${it.name}" }
_out3.println { "match: ${it.name}" }
