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

package nextflow.extension

import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
/**
 * Implements {@link DataflowExtensions#toList(groovyx.gpars.dataflow.DataflowReadChannel)}  operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ToListOp {

    static DataflowVariable apply( DataflowReadChannel source ) {
        assert source != null
        final target = new DataflowVariable()
        DataflowExtensions.reduceImpl(source, target, []) { list, item -> list << item }
        return target
    }

}
