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

package nextflow.script


import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.Global
import nextflow.Session

/**
 * Holds the current execution context
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ExecutionStack {

    static private List<ExecutionContext> stack = new ArrayList<>()

    static ExecutionContext current() {
        stack ? stack.get(0) : null
    }

    static boolean withinWorkflow() {
        for( def entry : stack ) {
            if( entry instanceof WorkflowDef )
                return true
        }
        return false
    }

    static WorkflowBinding binding() {
        stack ? current().getBinding() : (Global.session as Session).getBinding()
    }

    static BaseScript script() {
        for( def item in stack ) {
            if( item instanceof BaseScript )
                return item
        }
        throw new IllegalStateException("Missing execution script")
    }

    static BaseScript owner() {
        def c = current()
        if( c instanceof BaseScript )
            return c
        if( c instanceof WorkflowDef )
            return c.getOwner()
        throw new IllegalStateException("Not a valid scope object: [${c.getClass().getName()}] $this")
    }

    static WorkflowDef workflow() {
        final ctx = current()
        ctx instanceof WorkflowDef ? ctx : null
    }

    static void push(ExecutionContext script) {
        stack.push(script)
    }

    static ExecutionContext pop() {
        stack.pop()
    }

    static int size() {
        stack.size()
    }

    @PackageScope
    static void reset() {
        stack = new ArrayList<>()
    }

}
