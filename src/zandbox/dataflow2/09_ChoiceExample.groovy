/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
 * --== Choices ==--
 *
 *  The binaryChoice() and choice() methods allow you to send a value to one out of two (or many) output channels,
 *  as indicated by the return value from a closure
 */

queue1 = new DataflowQueue<>()
queue2 = new DataflowQueue<>()
queue3 = new DataflowQueue<>()
queue4 = new DataflowQueue<>()

queue1 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9
queue1.choice([queue2, queue3, queue4]) {a -> a % 3}
//queue1.binaryChoice(queue2, queue3) {a -> a > 0}

println "queue2: ${queue2.getVal()} - ${queue2.getVal()} - ${queue2.getVal()} "
println "queue3: ${queue3.getVal()} - ${queue3.getVal()} - ${queue3.getVal()} "
println "queue4: ${queue4.getVal()} - ${queue4.getVal()} - ${queue4.getVal()} "