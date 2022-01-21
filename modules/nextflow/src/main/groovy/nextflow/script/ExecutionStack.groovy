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

import java.nio.file.Path

/**
 * Holds the current execution context
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ExecutionStack {

    static private List<ExecutionContext> stack = new ArrayList<>()

    @TupleConstructor
    private static class TraceElement {
        int idx
        ComponentDef called
    }

    protected static final int CALL_TRACE_LIMIT = 40
    protected static int currentCallIdx = 1
    protected static LinkedList<TraceElement> callTrace = new LinkedList<TraceElement>()

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

    static int pushCallTrace(ComponentDef called) {
        callTrace.addLast(new TraceElement(currentCallIdx, called))
        if (callTrace.size() > CALL_TRACE_LIMIT)
            callTrace.removeFirst()
        return currentCallIdx++
    }

    static void popCallTrace(int idx, ComponentDef called) {
        callTrace.addLast(new TraceElement(-idx, called))
        if (callTrace.size() > CALL_TRACE_LIMIT)
            callTrace.removeFirst()
    }

    static Throwable registeredException

    static registerStackException(Throwable e) {
        if (e !== registeredException) {
            registeredException = e
            def last_cause = e
            while (last_cause.getCause() != null)
                last_cause = last_cause.getCause()
            last_cause.initCause(new ErrorInNextflowPipeline(e))
        }
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

@CompileStatic
class ErrorInNextflowPipeline extends Throwable {

    Throwable source;

    ErrorInNextflowPipeline(Throwable source) {
        this.source = source
        this.setStackTrace(new StackTraceElement[0])
    }

    static formatCalled(ComponentDef called) {
        "${called.name}(${called.getClass().getSimpleName()})"
    }

    static recFindUserStack(Map<String,Path> scriptPaths, List builtStack, Throwable e) {
        if (e.getCause())
            recFindUserStack(scriptPaths, builtStack, e.getCause())
        e.getStackTrace().each {
            def m = scriptPaths[it.className.replaceFirst(/\$.*$/,'')]
            if (m) {
                def lineNb = it.getLineNumber()
                builtStack.push("${it.getMethodName()}"
                        + "(${Global.session.baseDir.toAbsolutePath().relativize(m.toAbsolutePath())}"
                        + ":line ${lineNb > 0 ? lineNb : '?'})")
            }
        }
    }

    @Override
    String getMessage() {
        String res = "\n\nCALL HISTORY:\n" +
                (ExecutionStack.currentCallIdx > ExecutionStack.CALL_TRACE_LIMIT ? "\t..." : "")
        int last_idx = -1
        ExecutionStack.callTrace.each {
            if (it.idx < 0) {
                res += it.called == null ? "\n\t!${-it.idx} interrupted"
                        : it.idx == -last_idx ? '}' : "\n\t}${-it.idx}"
            } else {
                res += "\n\t#${it.idx}  ${formatCalled(it.called)} {"
            }
            last_idx = it.idx
        }
        res += "\n\nCALL STACK:\n"
        def stack = []
        recFindUserStack(ScriptMeta.allScriptNames(), stack, source)
        def last = stack.size()-1
        stack.eachWithIndex{ def line, int i ->
            res += '\t' + ('  ' * i) + "at ${line}" + (i == last ? '  <- HERE\n' : '\n')
        }
        res += "\n\n"
    }
}
