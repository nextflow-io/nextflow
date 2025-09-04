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
 */

package nextflow.script

import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.Nextflow

/**
 * Helper class for process entry execution feature.
 * 
 * This feature enables direct execution of Nextflow processes without explicit workflows:
 * - Single process scripts run automatically: `nextflow run script.nf --param value`
 * - Multi-process scripts use entry selection: `nextflow run script.nf -entry process:name --param value`
 * - Command-line parameters are mapped directly to process inputs
 * - Supports all standard Nextflow input types: val, path, env, tuple, each
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ProcessEntryHandler {

    private final BaseScript script
    private final Session session
    private final ScriptMeta meta

    ProcessEntryHandler(BaseScript script, Session session, ScriptMeta meta) {
        this.script = script
        this.session = session
        this.meta = meta
    }

    /**
     * Creates a workflow to execute a single standalone process automatically.
     * This allows single-process scripts to run without requiring the -entry option.
     * 
     * @return WorkflowDef that executes the single process with parameter mapping
     */
    WorkflowDef createSingleProcessWorkflow() {
        def processNames = meta.getLocalProcessNames()
        if( processNames.size() != 1 ) {
            throw new IllegalStateException("Expected exactly one process, found: ${processNames.size()}")
        }
        
        final processName = processNames.first()
        final processDef = meta.getProcess(processName)
        
        return createProcessWorkflow(processDef)
    }

    /**
     * Creates a workflow to execute a specific process with parameter mapping.
     * This enables process execution via the -entry process:NAME syntax.
     * 
     * @param processDef The ProcessDef object for the target process
     * @return WorkflowDef that executes the process with parameter mapping
     */
    WorkflowDef createProcessEntryWorkflow(ProcessDef processDef) {
        return createProcessWorkflow(processDef)
    }

    /**
     * Creates a workflow to execute the specified process with automatic parameter mapping.
     */
    private WorkflowDef createProcessWorkflow(ProcessDef processDef) {
        final processName = processDef.name
        
        // Create the workflow execution logic
        def workflowLogic = { ->
            // Get input parameter values and execute the process
            def inputArgs = getProcessInputArguments(processDef)
            script.invokeMethod(processName, inputArgs as Object[])
        }
        
        // Create workflow metadata
        def sourceCode = "    // Auto-generated process workflow\n    ${processName}(...)"
        
        // Wrap in BodyDef closure as expected by WorkflowDef constructor
        def workflowBody = { ->
            return new BodyDef(workflowLogic, sourceCode, 'workflow')
        }
        
        return new WorkflowDef(script, workflowBody)
    }

    /**
     * Gets the input arguments for a process by parsing input parameter names
     * and looking up corresponding values from session.params.
     * 
     * @param processDef The ProcessDef object containing the process definition
     * @return List of parameter values to pass to the process
     */
    private List getProcessInputArguments(ProcessDef processDef) {
        try {
            def inputNames = parseProcessInputNames(processDef)
            
            if( inputNames.isEmpty() ) {
                return []
            }
            
            // Map parameter names to values from session.params
            List inputArgs = []
            for( String paramName : inputNames ) {
                def paramValue = session.params.get(paramName)
                
                if( paramValue != null ) {
                    // Convert string paths to Path objects using file() helper
                    if( paramValue instanceof String && (paramValue.startsWith('/') || paramValue.contains('.'))) {
                        paramValue = Nextflow.file(paramValue)
                    }
                    inputArgs.add(paramValue)
                } else {
                    throw new IllegalArgumentException("Missing required parameter: --${paramName}")
                }
            }
            
            return inputArgs
            
        } catch (Exception e) {
            log.error "Failed to get input arguments for process ${processDef.name}: ${e.message}"
            throw e
        }
    }
    
    /**
     * Parses the process body to extract input parameter names by intercepting
     * Nextflow's internal compiled method calls (_in_val, _in_path, etc.).
     *
     * @param processDef The ProcessDef containing the raw process body
     * @return List of input parameter names found in the process
     */
    private List<String> parseProcessInputNames(ProcessDef processDef) {
        def inputNames = []
        
        // Create delegate to capture Nextflow's internal input method calls
        def delegate = new Object() {
            def _in_val(tokenVar) { inputNames.add(tokenVar.name.toString()) }
            def _in_path(tokenVar) { inputNames.add(tokenVar.name.toString()) }
            def _in_file(tokenVar) { inputNames.add(tokenVar.name.toString()) }
            def _in_env(tokenVar) { inputNames.add(tokenVar.name.toString()) }
            def _in_each(tokenVar) { inputNames.add(tokenVar.name.toString()) }
            
            def _in_tuple(Object... items) {
                for( item in items ) {
                    if( item?.hasProperty('name') ) {
                        inputNames.add(item.name.toString())
                    }
                }
            }
            
            // Handle legacy input block syntax for backward compatibility
            def input(Closure inputBody) {
                def inputDelegate = new Object() {
                    def val(name) { inputNames.add(name.toString()) }
                    def path(name) { inputNames.add(name.toString()) }
                    def file(name) { inputNames.add(name.toString()) }
                    def env(name) { inputNames.add(name.toString()) }
                    def each(name) { inputNames.add(name.toString()) }
                    def tuple(Object... items) {
                        for( item in items ) {
                            if( item instanceof String || (item instanceof groovy.lang.GString) ) {
                                inputNames.add(item.toString())
                            }
                        }
                    }
                    def methodMissing(String name, args) {
                        for( arg in args ) {
                            if( arg instanceof String || (arg instanceof groovy.lang.GString) ) {
                                inputNames.add(arg.toString())
                            }
                        }
                    }
                }
                inputBody.delegate = inputDelegate
                inputBody.resolveStrategy = Closure.DELEGATE_FIRST
                inputBody.call()
            }
            
            // Ignore all other method calls during parsing
            def methodMissing(String name, args) { /* ignore */ }
        }
        
        // Execute the process body with our capturing delegate
        def bodyClone = processDef.rawBody.clone()
        bodyClone.delegate = delegate
        bodyClone.resolveStrategy = Closure.DELEGATE_FIRST
        
        try {
            bodyClone.call()
        } catch (Exception e) {
            // Ignore exceptions during parsing - we only want to capture input names
        }
        
        return inputNames
    }
}