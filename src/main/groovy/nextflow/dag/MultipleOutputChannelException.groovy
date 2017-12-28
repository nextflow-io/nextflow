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

package nextflow.dag

import groovy.transform.PackageScope

/**
 * Exception raised then the same channel is declared as output more
 * than one time
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

class MultipleOutputChannelException extends Exception {

    MultipleOutputChannelException( String name, DAG.Vertex duplicate, DAG.Vertex existing ) {
        super(message(name,duplicate,existing))
    }

    @PackageScope
    static message( String name, DAG.Vertex duplicate, DAG.Vertex existing ) {
        if( !name ) {
            return 'Channels cannot be used as output in more than one process or operator'
        }

        if( !duplicate || duplicate.type != DAG.Type.PROCESS ) {
            return "Channel `$name` has been used as an output by more than a process or an operator"
        }

        String message = "Channel `$name` has been used twice as an output by process `${duplicate.label}`"
        if( existing != duplicate )  {
            message += existing?.type == DAG.Type.PROCESS  ? " and process `${existing.label}`" : " and another operator"
        }

        return message
    }
}
