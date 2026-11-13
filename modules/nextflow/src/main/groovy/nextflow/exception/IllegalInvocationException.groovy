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
