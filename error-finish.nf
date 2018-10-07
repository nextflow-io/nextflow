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


process foo {
  echo true
  errorStrategy 'finish'
  input:  each x from 1,2,3
  output: stdout into results

  script:
  if( x != 3 )
  """
    echo run_$x  
    sleep 5
  """
  else
  """
    exit 99
  """
}

process bar {
  input:  file 'x' from results

  script:
  '''
  cat x
  '''
}


workflow.onError {
  println "success: $workflow.success"
  println "exitStatus: $workflow.exitStatus"
}