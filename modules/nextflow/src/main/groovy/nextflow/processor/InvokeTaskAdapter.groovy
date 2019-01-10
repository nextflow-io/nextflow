/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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