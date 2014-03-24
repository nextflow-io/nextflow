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
 * --==  Tapping into the pipeline ==--
 *
 * + Like split() the tap() method allows you to fork the data flow into multiple channels.
 *   Tapping, however, is slightly more convenient in some scenarios, since it treats one of
 *   the two new forks as the successor of the pipeline.
 */

def queue = new DataflowQueue()
def logChannel = new DataflowQueue()
def printChannel = new DataflowQueue()

queue.chainWith {it * 2}.chainWith{it + 10}.tap(logChannel).into(printChannel)

queue << 1 << 2 << 3

logChannel.each { println "log: $it"  }

//printChannel.each { println "print: $it"  }

sleep 1000