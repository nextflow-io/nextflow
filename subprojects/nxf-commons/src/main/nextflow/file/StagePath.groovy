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

package nextflow.file

import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.Paths

import org.codehaus.groovy.runtime.InvokerHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Deprecated
class StagePath {

    final Path target

    StagePath( Path source ) {
        this.target = source
    }

    StagePath( String fileName ) {
        this.target = Paths.get(fileName)
    }

    Object getProperty( String name ) {
        InvokerHelper.getProperty(target,name)
    }

    void setProperty(String property, Object newValue) {
        InvokerHelper.setProperty(target, property, newValue)
    }

    Object invokeMethod(String name, Object args) {
        InvokerHelper.invokeMethod(target, name, args)
    }

    Path getFileName() {
        target.getFileName()
    }

    Path toRealPath(LinkOption... options) {
        return target
    }

    String toString() {
        target.getFileName().toString()
    }
}
