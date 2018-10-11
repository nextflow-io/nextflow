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

package nextflow.processor

import groovy.transform.CompileStatic

/**
 * Adapter closure to call the {@link TaskProcessor#invokeTask(java.lang.Object)} method
 */
@CompileStatic
class InvokeTaskAdapter extends Closure {

    private int numOfParams

    private TaskProcessor processor

    InvokeTaskAdapter(TaskProcessor p, int n) {
        super(null, null);
        processor = p
        numOfParams = n
    }

    @Override
    int getMaximumNumberOfParameters() {
        numOfParams
    }

    @Override
    Class[] getParameterTypes() {
        def result = new Class[numOfParams]
        for( int i=0; i<result.size(); i++ )
            result[i] = Object
        return result
    }

    @Override
    Object call(final Object arguments) {
        processor.invokeTask(arguments)
        return null
    }

    @Override
    Object call(final Object... args) {
        processor.invokeTask(args as List)
        return null
    }

    @Override
    Object call() {
        throw new UnsupportedOperationException()
    }
}