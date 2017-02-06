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

import static groovyx.gpars.dataflow.Dataflow.task

import groovyx.gpars.dataflow.DataflowBroadcast

def broadcastStream = new DataflowBroadcast()
def stream1 = broadcastStream.createReadChannel()
def stream2 = broadcastStream.createReadChannel()


task {

    while( true ) {
        println "Channel1: " + stream1.val
    }

}

task {

    while( true ) {
        println "Channel2: " + stream2.val
    }

}



broadcastStream << 'Message1'
broadcastStream << 'Message2'
broadcastStream << 'Message3'

// sleep otherwise it terminates before printout the result
sleep 1000

