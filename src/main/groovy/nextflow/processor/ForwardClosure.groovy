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
import groovyx.gpars.dataflow.operator.DataflowProcessor

/**
 * Implements the closure which *combines* all the iteration
 *
 * @param numOfInputs Number of in/out channel
 * @param indexes The list of indexes which identify the position of iterators in the input channels
 * @return The closure implementing the iteration/forwarding logic
 */
@CompileStatic
class ForwardClosure extends Closure {

    final private Integer len

    final private int numOfParams

    final private List<Integer> indexes

    ForwardClosure(int len, List<Integer> indexes) {
        super(null, null);
        this.len = len
        this.numOfParams = len+1
        this.indexes = indexes
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
    Object call(final Object... args) {
        final target = ((DataflowProcessor) getDelegate())
        /*
         * Explaining the following code:
         *
         * - 'out' holds the list of all input values which need to be forwarded (bound) to output as many
         *   times are the items in the iteration list
         *
         * - the iteration list(s) is (are) passed like in the closure inputs like the other values,
         *   the *indexes* argument defines the which of them are the iteration lists
         *
         * - 'itr' holds the list of all iteration lists
         *
         * - using the groovy method a combination of all values is create (cartesian product)
         *   see http://groovy.codehaus.org/groovy-jdk/java/util/Collection.html#combinations()
         *
         * - the resulting values are replaced in the 'out' array of values and forwarded out
         *
         */

        def out = args[0..-2]
        def itr = indexes.collect { args[it] }
        List<List> cmb = itr.combinations()

        for( int i=0; i<cmb.size(); i++ ) {
            List entries = cmb[i]
            int count = 0
            for( int j=0; j<len; j++ ) {
                if( j in this.indexes ) {
                    out[j] = entries[count++]
                }
            }

            target.bindAllOutputValues( out as Object[] )
        }

        return null
    }

    @Override
    Object call(final Object arguments) {
        throw new UnsupportedOperationException()
    }

    @Override
    Object call() {
        throw new UnsupportedOperationException()
    }
}