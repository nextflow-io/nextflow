/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.extension.op

import nextflow.prov.OperatorRun

/**
 * A closure that wraps the execution of an operator target code (closure)
 * and maps the inputs and outputs to the corresponding operator run.
 *
 * This class extends {@link OpClosure} assuming that all results are `"accumulated"`
 * as they were executed in the same operator run
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class OpAccumulatorClosure extends OpClosure {

    OpAccumulatorClosure(Closure code) {
        super(code)
        setPreviousRun(new OperatorRun())
    }

    @Override
    protected OperatorRun runInstance() {
        return getPreviousRun()
    }
}
