// GPars - Groovy Parallel Systems
//
// Copyright © 2008-11  The original author or authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Modified by Seqera Labs for Nextflow (2026): updated the class description; moved from groovyx.gpars.dataflow.operator to Nextflow, its only user.

package nextflow.extension;

import groovy.lang.Closure;
import groovyx.gpars.dataflow.operator.DataflowProcessor;
import org.codehaus.groovy.runtime.NullObject;

/**
 * Wraps the supplied closure inside the body of an operator.
 * This wrapper closure is responsible for writing the calculated result into the operator's output through its bind() methods.
 *
 * @author Vaclav Pech
 */

public final class ChainWithClosure<V> extends Closure {
    private final Closure code;

    public ChainWithClosure(final Closure code) {
        super(null, null);
        this.code = code;
    }

    @Override
    public int getMaximumNumberOfParameters() {
        return code.getMaximumNumberOfParameters();
    }

    @Override
    public Class[] getParameterTypes() {
        return code.getParameterTypes();
    }

    @Override
    public void setDelegate(final Object delegate) {
        super.setDelegate(delegate);
        code.setDelegate(delegate);
    }

    @Override
    public void setResolveStrategy(final int resolveStrategy) {
        super.setResolveStrategy(resolveStrategy);
        code.setResolveStrategy(resolveStrategy);
    }

    @Override
    public Object call(final Object arguments) {
        final V result = (V) code.call(arguments);
        if (result != NullObject.getNullObject()) ((DataflowProcessor) getDelegate()).bindAllOutputsAtomically(result);
        return result;
    }

    @Override
    public Object call(final Object... args) {
        final V result = (V) code.call(args);
        if (result != NullObject.getNullObject()) ((DataflowProcessor) getDelegate()).bindAllOutputsAtomically(result);
        return result;
    }

    @Override
    public Object call() {
        final V result = (V) code.call();
        if (result != NullObject.getNullObject()) ((DataflowProcessor) getDelegate()).bindAllOutputsAtomically(result);
        return result;
    }
}
