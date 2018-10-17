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

x = 100
y = 200

process foo {
  input:
  file fastq from 'dummy'

  output:
  val 'Hello' into str_channel
  val "${fastq.baseName}-${x}.out" into exp_channel
  val x into x_channel
  val y into y_channel

  script:
  y = 'two hundred'
  """
  echo bar
  """
}

x_channel.println { "x: $it" }
y_channel.println { "y: $it" }
str_channel.println { "str: $it" }
exp_channel.println { "exp: $it" }
