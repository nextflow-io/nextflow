/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.script.params

import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.InheritConstructors
import groovy.transform.Memoized

/**
 * Model process `output: eval PARAM` definition
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class CmdEvalParam extends BaseOutParam implements OptionalParam {

    private static AtomicInteger counter = new AtomicInteger()

    private Object target

    private int count

    {
        count = counter.incrementAndGet()
    }

    String getName() {
        return "nxf_out_eval_${count}"
    }

    BaseOutParam bind( def obj ) {
        if( obj !instanceof CharSequence )
            throw new IllegalArgumentException("Invalid argument for command output: $this")
        // the target value object
        target = obj
        return this
    }

    @Memoized
    String getTarget(Map<String,Object> context) {
        return target instanceof GString
            ? target.cloneAsLazy(context).toString()
            : target.toString()
    }
}
