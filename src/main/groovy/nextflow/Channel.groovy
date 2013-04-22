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

package nextflow
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.operator.PoisonPill

/**
 * A DataflowQueue with that may be pre-loaded with a finite number of items
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class Channel<T> extends DataflowQueue<T> {

    Channel(Collection<T> values) {
        if ( values )  {
            values.each { this << it }
            // since queue is 'finite' close it by a poison pill
            // so the operator will stop on when all values in the queue are consumed
            // (otherwise it will wait forever for a new entry)
            this << PoisonPill.instance
        }
    }

    Channel( T ... values ) {
        this( values ? values as List : [] )
    }
}
