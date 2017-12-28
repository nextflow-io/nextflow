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

package nextflow.extension
import groovyx.gpars.dataflow.operator.DataflowProcessor
/**
 * Implement the default `merge` operator logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultMergeClosure extends Closure {

    private int numOfParams

    DefaultMergeClosure(int n) {
        super(null, null);
        numOfParams = n
    }

    @Override
    public int getMaximumNumberOfParameters() {
         numOfParams
    }

    @Override
    public Class[] getParameterTypes() {
        Collections.nCopies(numOfParams, Object)
    }

    @Override
    public void setDelegate(final Object delegate) {
        super.setDelegate(delegate);
    }

    @Override
    public void setResolveStrategy(final int resolveStrategy) {
        super.setResolveStrategy(resolveStrategy);
    }

    @Override
    public Object call(final Object arguments) {
        throw new UnsupportedOperationException()
    }

    @Override
    public Object call(final Object... args) {
        final result = []
        for( int i=0; i<args.size(); i++ )
            DataflowHelper.addToList(result, args[i])
        ((DataflowProcessor) getDelegate()).bindAllOutputsAtomically(result);
        return result;
    }

    @Override
    public Object call() {
        throw new UnsupportedOperationException()
    }

}
