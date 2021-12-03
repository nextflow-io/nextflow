/*
 * Copyright 2020-2021, Seqera Labs
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
import groovy.transform.TupleConstructor
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

    protected static List<ComponentDef> fullCallStack = new ArrayList<ComponentDef>()

    @TupleConstructor
    private static class TraceElement {
        int depth
        ComponentDef called
    }
    private static LinkedList<TraceElement> callTrace = new LinkedList<TraceElement>()

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

    static List<WorkflowDef> workflowStack() {
        return (List<WorkflowDef>)(Object)stack.findAll {it instanceof WorkflowDef}
    }

    static WorkflowDef workflow() {
        final ctx = current()
        ctx instanceof WorkflowDef ? ctx : null
    }

    static void pushFull(ComponentDef called) {
        callTrace.addLast(new TraceElement(fullCallStack.size(), called))
        fullCallStack.add(called)
    }

    static void popFull() {
        fullCallStack.pop()
    }

    static List<ComponentDef> registeredStack
    static Throwable registeredException

    static registerStackException(Throwable e) {
        if (e !== registeredException) {
            registeredException = e
            e.addSuppressed(new NEXTFLOW_PIPELINE_STACK())
            registeredStack = (List<ComponentDef>)fullCallStack.clone()
        }
    }

    static void push(ExecutionContext script) {
        if (script instanceof ComponentDef)
            pushFull(script)
        stack.push(script)
    }

    static ExecutionContext pop() {
        def res = stack.pop()
        if (res instanceof ComponentDef)
            fullCallStack.pop()
        res
    }

    static int size() {
        stack.size()
    }

    @PackageScope
    static void reset() {
        stack = new ArrayList<>()
    }

}

@CompileStatic
class NEXTFLOW_PIPELINE_STACK extends Throwable {

    NEXTFLOW_PIPELINE_STACK() {
        this.setStackTrace(new StackTraceElement[0])
    }
    @Override
    String getMessage() {
        String res = "\n"
        ExecutionStack.fullCallStack.eachWithIndex {called, idx ->
            res += "${"\t"*idx}${called.name}(${called.getClass().toString().replace("Def","")})\n"
        }
        res
    }
}
