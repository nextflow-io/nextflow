/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

package nextflow.dag
import groovy.transform.ToString
import nextflow.dag.DAG.ChannelHandler
import nextflow.dag.DAG.Vertex
/**
 * Exception thrown when the same channel is declared as input more than one time.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true)
class MultipleInputChannelException extends Exception {

    String name

    ChannelHandler channel

    Vertex duplicate

    Vertex existing

    MultipleInputChannelException( String name, ChannelHandler channel, Vertex duplicate, Vertex existing ) {
        super(message(name,duplicate,existing))
        this.name = name
        this.channel = channel
        this.duplicate = duplicate
        this.existing = existing
    }


    private static String message(String name, Vertex duplicate, Vertex existing) {
        if( !name ) {
            return 'Channels cannot be used as input in more than one process or operator'
        }

        if( duplicate.type != DAG.Type.PROCESS ) {
            return "Channel `$name` has been used as an input by more than a process or an operator"
        }

        String message = "Channel `$name` has been used twice as an input by process `${duplicate.label}`"
        if( existing.type == DAG.Type.PROCESS )
            message += " and process `${existing.label}`"
        else
            message += " and another operator"
        return message
    }

}
