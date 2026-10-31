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

/**
 * Identity closure used with {@link ChainWithClosure} by the copying operators
 * to forward each input value unchanged to the operator's outputs.
 *
 * @author Vaclav Pech
 */

public final class CopyChannelsClosure<V> extends Closure {

    private static final Class[] PARAMETER_TYPES = {Object.class};

    public CopyChannelsClosure() {
        super(null, null);
    }

    @Override
    public int getMaximumNumberOfParameters() {
        return 1;
    }

    @SuppressWarnings({"ReturnOfCollectionOrArrayField"})
    @Override
    public Class[] getParameterTypes() {
        return PARAMETER_TYPES;
    }

    @Override
    public Object call(final Object arguments) {
        return arguments;
    }

    @SuppressWarnings({"OverloadedVarargsMethod"})
    @Override
    public Object call(final Object... args) {
        return args[0];
    }

    @Override
    public Object call() {
        return null;
    }
}
