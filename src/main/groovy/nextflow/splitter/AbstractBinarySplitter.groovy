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

package nextflow.splitter

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
/**
 * Base class for splitter handling binary data
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
@InheritConstructors
abstract class AbstractBinarySplitter extends AbstractSplitter<InputStream> {

    protected InputStream normalizeSource( obj ) {

        if( obj instanceof InputStream )
            return (InputStream) obj

        if( obj instanceof byte[] )
            return new ByteArrayInputStream((byte[])obj)

        if( obj instanceof CharSequence )
            return new ByteArrayInputStream(obj.toString().bytes)

        if( obj instanceof Path )
            return newInputStream(obj)

        if( obj instanceof File )
            newInputStream(obj.toPath())

        throw new IllegalAccessException("Object of class '${obj.class.name}' does not support 'splitter' methods")

    }



}
