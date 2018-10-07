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


foo = Channel.from('abc')
loop = Channel.create()

process bar {
  input:
  file x from foo .mix(loop.take(2))
  output:
  file q

  """
  rev $x > q
  """
}

process baz {
  input:
  file q from q
  output:
  file z into loop
  file z into result

  """
  cat $q > z
  echo 'hello' >> z
  """
}

result.last().println { it.text }
