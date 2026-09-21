/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.script.params

import groovy.transform.InheritConstructors
import groovy.transform.Memoized

/**
 * Model process `output: eval PARAM` definition
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class CmdEvalParam extends BaseOutParam implements OptionalParam {

    private Object target

    /**
     * The eval output name is derived from the parameter position in the process
     * output block, so that it is stable across runs. Note: it must not depend on
     * a global counter, otherwise adding or removing an `eval` output in one process
     * would change the name -- and therefore the task hash -- of unrelated processes,
     * invalidating their cache. See https://github.com/nextflow-io/nextflow/issues/7572
     */
    String getName() {
        return mapIndex >= 0
            ? "nxf_out_eval_${index}_${mapIndex}"
            : "nxf_out_eval_${index}"
    }

    BaseOutParam bind( def obj ) {
        if( obj !instanceof CharSequence && obj !instanceof Closure )
            throw new IllegalArgumentException("Invalid argument for command output: $this")
        // the target value object
        target = obj
        return this
    }

    @Memoized
    String getTarget(Map<String,Object> context) {
        return target instanceof Closure
            ? context.with(target)
            : target instanceof GString
            ? target.cloneAsLazy(context).toString()
            : target.toString()
    }
}
