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

package nextflow.file

import java.nio.file.Path

import groovy.transform.CompileStatic

/**
 * A Path object whose {@link Object#toString()} method return *always* the relative file name.
 * <p>
 *   This is needed when interpolating file in the process scripts
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class StagePath implements Path {

    @Delegate
    Path target

    StagePath(Path path) {
        assert path
        this.target = path
    }

    StagePath(String name) {
        this.target = FileHelper.asPath(name)
    }

    String toString() {
        return target.getFileName().toString()
    }

}
