/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

import groovyx.gpars.dataflow.DataflowQueue

/**
 * --== Filtering ==-
 *
 * + The filter() method allows to filter data in the pipeline using boolean predicates.
 *
 */


queue1 = new DataflowQueue()
queue2 = queue1.filter {num -> num % 2 != 0 }

// alternative version:
// queue1.filter(odd) into ( queue2 = new DataflowQueue() )

(1..5).each {queue1 << it}

assert 1 == queue2.val
assert 3 == queue2.val
assert 5 == queue2.val