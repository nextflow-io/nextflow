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

package nextflow.util

import java.lang.reflect.Field

import groovy.transform.CompileStatic

/**
 * Hack the Java process to be able to access the Unix process id
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class PosixProcess {

    @Delegate
    Process target

    PosixProcess( Process process ) {
        assert process

        target = process
    }

    /**
     * Access the *target* Unix process pid, via a reflection trick
     */
    @Lazy Integer pid = {

        Field field = target.class.getDeclaredField("pid");
        field.setAccessible(true);
        return field.getInt(target);

    } ()


    /**
     * Check if the submitted job has terminated its execution
     */
    boolean hasExited() {

        Field field = target.class.getDeclaredField("hasExited");
        field.setAccessible(true);
        return field.getBoolean(target);

    }

}
