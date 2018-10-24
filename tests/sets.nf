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

x = Channel.from( ['a', 'file1'], ['b','file2'] )

process touch {

  input:
    set ( id, fileName ) from x
  output:
    set ( id, 'file*' ) into z


  /
  echo Creating $id
  touch $fileName
  /
}

process makeFiles {
  input:
    set( id, 'file_x' ) from z

  output:
    set( id, '*') into q mode flatten

  /
   cp file_x copy_$id
   touch beta_$id
  /

}

q.subscribe { println it }
