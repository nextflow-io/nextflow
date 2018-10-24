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
  * This test verify the `workflow.onComplete` defined in the config file
  * works OK.
  *
  * See file `$PWD/checks/workflow-oncomplete.nf/.config` for details
  */

echo true

params.command = 'echo'
cheers = Channel.from 'Bojour', 'Ciao', 'Hello', 'Hola', 'Γεια σου'

process sayHello {
  input: 
  val x from cheers
  
  """
  ${params.command} '$x world!'
  """
}