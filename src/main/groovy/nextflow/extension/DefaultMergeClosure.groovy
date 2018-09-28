/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
