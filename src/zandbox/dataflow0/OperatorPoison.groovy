/*
 * Copyright (c) 2012, the authors.
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

import static groovyx.gpars.dataflow.Dataflow.operator

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
/**
 * Shows how to build operators using the ProcessingNode class
 */

final DataflowQueue aValues = new DataflowQueue()

//Now the operator is running and processing the data
aValues << 10 << 20 << 30 << 40 << 50 << PoisonPill.instance


def listener = new DataflowEventAdapter() {
    @Override
    public boolean onException(final DataflowProcessor processor, final Throwable e) {
        e.printStackTrace()
        return true;
    }

}


operator(inputs: [aValues], outputs:[], listeners: [listener])  { x ->
    println x
}.join()
