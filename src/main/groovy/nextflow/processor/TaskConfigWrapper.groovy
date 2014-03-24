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

package nextflow.processor

/**
 * Wrap a {@code TaskConfig} object in order to modify the behaviour of
 * {@code #getProperty} so that it raises a {@code MissingPropertyException}
 * when a property is not defined in the delegate map
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

class TaskConfigWrapper extends TaskConfig {

    TaskConfigWrapper( TaskConfig source ) {
        super(source)
    }

    def getProperty( String name ) {
        if( containsKey(name) ) {
            return super.getProperty(name)
        }

        throw new MissingPropertyException("Unknown variable '$name'")
    }

}
