package nextflow.exception

import nextflow.script.ComponentDef
import nextflow.script.ProcessDef
import nextflow.script.WorkflowDef

/**
 * Exception thrown when a module component is invoked
 * in a wrong context eg when invoking a process outside a workflow scope
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IllegalInvocationException extends ProcessException {

    IllegalInvocationException(ComponentDef component) {
        super(message(component))
    }

    static private String message(ComponentDef component) {
        if( component instanceof WorkflowDef )
            return "Workflow '$component.name' can only be invoked from workflow context"

        if( component instanceof ProcessDef )
            return "Process '$component.name' can only be invoked from a workflow context"

        return "Invalid invocation context: $component.name"
    }
}
