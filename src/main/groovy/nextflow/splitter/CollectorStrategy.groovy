/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

package nextflow.splitter

import groovy.transform.CompileStatic

/**
 * Defines an abstract contract to collect items produced by splitters into chunks
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
interface CollectorStrategy {

    /**
     * @param item Add an item to the current chunk
     */
    void add( Object item )

    /**
     * Skip to the next chunk. All following {@link #add(java.lang.Object)} invocations will
     * refer to a new chunk
     */
    void next()

    /**
     * @return {@code true} if at least an item as been added
     */
    boolean hasChunk()

    /**
     * @return The chunk made up of the added items
     */
    Object getChunk()

}
