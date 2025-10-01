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
 * - Multi-process scripts run the first process automatically: `nextflow run script.nf --param value`
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
     * Creates a workflow to execute a standalone process automatically.
     * For single process scripts, executes that process.
     * For multi-process scripts, executes the first process.
     * 
     * @return WorkflowDef that executes the process with parameter mapping
     */
    WorkflowDef createAutoProcessWorkflow() {
        def processNames = meta.getLocalProcessNames()
        if( processNames.isEmpty() ) {
            throw new IllegalStateException("No processes found for auto-execution")
        }
        
        // Always pick the first process (whether single or multiple processes)
        final processName = processNames.first()
        final processDef = meta.getProcess(processName)
        
        return createProcessWorkflow(processDef)
    }

    /**
     * Creates a workflow to execute the specified process with automatic parameter mapping.
     */
    private WorkflowDef createProcessWorkflow(ProcessDef processDef) {
        final processName = processDef.name
        
        // Create a simple workflow body that just executes the process
        def workflowBodyClosure = { ->
            // Create the workflow execution logic
            def workflowExecutionClosure = { ->
                // Get input parameter values and execute the process
                def inputArgs = getProcessInputArguments(processDef)
                def processResult = script.invokeMethod(processName, inputArgs as Object[])
                
                return processResult
            }
            
            // Create the body definition with execution logic
            def sourceCode = "    // Auto-generated process workflow\n    ${processName}(...)"
            return new BodyDef(workflowExecutionClosure, sourceCode, 'workflow')
        }
        
        // Create a simple workflow definition without complex emission logic
        return new WorkflowDef(script, workflowBodyClosure)
    }

    /**
     * Gets the input arguments for a process by parsing input parameter structures
     * and mapping them from session.params, supporting dot notation for complex inputs.
     * 
     * @param processDef The ProcessDef object containing the process definition
     * @return List of parameter values to pass to the process
     */
    private List getProcessInputArguments(ProcessDef processDef) {
        try {
            log.debug "Getting input arguments for process: ${processDef.name}"
            log.debug "Session params: ${session.params}"
            
            def inputStructures = parseProcessInputStructures(processDef)
            log.debug "Parsed input structures: ${inputStructures}"
            
            if( inputStructures.isEmpty() ) {
                log.debug "No input structures found, returning empty list"
                return []
            }
            
            // Parse complex parameters from session.params (handles dot notation)
            def complexParams = parseComplexParameters(session.params)
            log.debug "Complex parameters: ${complexParams}"
            
            // Map input structures to actual values
            List inputArgs = []
            for( def inputDef : inputStructures ) {
                log.debug "Processing input definition: ${inputDef}"
                if( inputDef.type == 'tuple' ) {
                    // Handle tuple inputs - construct list with proper elements
                    List tupleElements = []
                    for( def element : inputDef.elements ) {
                        log.debug "Getting value for tuple element: ${element}"
                        def value = getValueForInput(element, complexParams)
                        tupleElements.add(value)
                    }
                    log.debug "Constructed tuple: ${tupleElements}"
                    inputArgs.add(tupleElements)
                } else {
                    // Handle simple inputs
                    def value = getValueForInput(inputDef, complexParams)
                    log.debug "Got simple input value: ${value}"
                    inputArgs.add(value)
                }
            }
            
            log.debug "Final input arguments: ${inputArgs}"
            return inputArgs
            
        } catch (Exception e) {
            log.error "Failed to get input arguments for process ${processDef.name}: ${e.message}"
            throw e
        }
    }
    
    /**
     * Parses the process body to extract input parameter structures by intercepting
     * Nextflow's internal compiled method calls (_in_val, _in_path, _in_tuple, etc.).
     *
     * @param processDef The ProcessDef containing the raw process body
     * @return List of input structures with type and name information
     */
    private List parseProcessInputStructures(ProcessDef processDef) {
        def inputStructures = []
        
        // Create delegate to capture Nextflow's internal input method calls
        def delegate = new Object() {
            def _in_val(tokenVar) { 
                def varName = extractVariableName(tokenVar)
                if( varName ) inputStructures.add([type: 'val', name: varName]) 
            }
            def _in_path(tokenVar) { 
                def varName = extractVariableName(tokenVar)
                if( varName ) inputStructures.add([type: 'path', name: varName]) 
            }
            def _in_file(tokenVar) { 
                def varName = extractVariableName(tokenVar)
                if( varName ) inputStructures.add([type: 'file', name: varName]) 
            }
            def _in_env(tokenVar) { 
                def varName = extractVariableName(tokenVar)
                if( varName ) inputStructures.add([type: 'env', name: varName]) 
            }
            def _in_each(tokenVar) { 
                def varName = extractVariableName(tokenVar)
                if( varName ) inputStructures.add([type: 'each', name: varName]) 
            }
            
            def extractVariableName(token) {
                if( token?.hasProperty('name') ) {
                    return token.name.toString()
                } else {
                    // Try to extract from string representation
                    def match = token.toString() =~ /TokenVar\(([^)]+)\)/
                    return match ? match[0][1] : null
                }
            }
            
            def _in_tuple(Object... items) {
                def tupleElements = []
                for( item in items ) {
                    log.debug "Processing tuple item: ${item} of class ${item?.getClass()?.getSimpleName()}"
                    
                    def itemType = 'val' // default
                    def itemName = null
                    
                    // Handle different token call types by checking class name
                    def className = item.getClass().getSimpleName()
                    if( className == 'TokenValCall' ) {
                        itemType = 'val'
                        itemName = extractVariableNameFromToken(item)
                    } else if( className == 'TokenPathCall' || className == 'TokenFileCall' ) {
                        itemType = 'path'
                        itemName = extractVariableNameFromToken(item)
                    } else if( className == 'TokenEnvCall' ) {
                        itemType = 'env'
                        itemName = extractVariableNameFromToken(item)
                    } else if( className == 'TokenEachCall' ) {
                        itemType = 'each'
                        itemName = extractVariableNameFromToken(item)
                    } else {
                        // Fallback: try to extract from string representation
                        if( item.toString().contains('TokenValCall') ) {
                            itemType = 'val'
                            def tokenVar = item.toString().find(/TokenVar\(([^)]+)\)/) { match, varName -> varName }
                            itemName = tokenVar
                        }
                    }
                    
                    if( itemName ) {
                        log.debug "Parsed tuple element: ${itemName} (${itemType})"
                        tupleElements.add([type: itemType, name: itemName])
                    } else {
                        log.warn "Could not parse tuple element: ${item} of class ${className}"
                    }
                }
                log.debug "Parsed tuple with ${tupleElements.size()} elements: ${tupleElements}"
                inputStructures.add([type: 'tuple', elements: tupleElements])
            }
            
            def extractVariableNameFromToken(token) {
                // Try to access the variable property directly
                try {
                    if( token.hasProperty('variable') && token.variable?.hasProperty('name') ) {
                        return token.variable.name.toString()
                    }
                    if( token.hasProperty('target') && token.target?.hasProperty('name') ) {
                        return token.target.name.toString()
                    }
                    if( token.hasProperty('name') ) {
                        return token.name.toString()
                    }
                    // Fallback to string parsing
                    def match = token.toString() =~ /TokenVar\(([^)]+)\)/
                    return match ? match[0][1] : null
                } catch( Exception e ) {
                    log.debug "Error extracting variable name from ${token}: ${e.message}"
                    return null
                }
            }
            
            // Handle legacy input block syntax for backward compatibility
            def input(Closure inputBody) {
                def inputDelegate = new Object() {
                    def val(name) { 
                        inputStructures.add([type: 'val', name: name.toString()]) 
                    }
                    def path(name) { 
                        inputStructures.add([type: 'path', name: name.toString()]) 
                    }
                    def file(name) { 
                        inputStructures.add([type: 'file', name: name.toString()]) 
                    }
                    def env(name) { 
                        inputStructures.add([type: 'env', name: name.toString()]) 
                    }
                    def each(name) { 
                        inputStructures.add([type: 'each', name: name.toString()]) 
                    }
                    def tuple(Object... items) {
                        def tupleElements = []
                        for( item in items ) {
                            if( item instanceof String || (item instanceof groovy.lang.GString) ) {
                                tupleElements.add([type: 'val', name: item.toString()])
                            }
                        }
                        inputStructures.add([type: 'tuple', elements: tupleElements])
                    }
                    def methodMissing(String name, args) {
                        for( arg in args ) {
                            if( arg instanceof String || (arg instanceof groovy.lang.GString) ) {
                                inputStructures.add([type: name, name: arg.toString()])
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
            // Ignore exceptions during parsing - we only want to capture input structures
        }
        
        return inputStructures
    }
    
    /**
     * Parses complex parameters with dot notation support.
     * Converts flat parameters like --meta.id=1 --meta.name=test to nested maps.
     * 
     * @param params Flat parameter map from session.params
     * @return Map with nested structures for complex parameters
     */
    private Map parseComplexParameters(Map params) {
        Map complexParams = [:]
        
        params.each { key, value ->
            def parts = key.toString().split('\\.')
            if( parts.length > 1 ) {
                // Handle dot notation - build nested map
                def current = complexParams
                for( int i = 0; i < parts.length - 1; i++ ) {
                    if( !current.containsKey(parts[i]) ) {
                        current[parts[i]] = [:]
                    }
                    current = current[parts[i]]
                }
                current[parts[-1]] = value
            } else {
                // Simple parameter
                complexParams[key] = value
            }
        }
        
        return complexParams
    }
    
    /**
     * Parses file input handling comma-separated values.
     * If the input contains commas, splits and returns a list of files.
     * Otherwise returns a single file object.
     * 
     * @param fileInput String representation of file path(s)
     * @return Single file object or List of file objects
     */
    protected Object parseFileInput(String fileInput) {
        if (fileInput.contains(',')) {
            // Split by comma, trim whitespace, and convert each to a file
            return fileInput.tokenize(',')
                .collect { it.trim() }
                .findAll { !it.isEmpty() }
                .collect { Nextflow.file(it) }
        } else {
            // Single file case - existing behavior
            return Nextflow.file(fileInput)
        }
    }
    
    /**
     * Gets the appropriate value for an input definition, handling type conversion.
     * 
     * @param inputDef Input definition with type and name
     * @param complexParams Parsed parameter map with nested structures
     * @return Properly typed value for the input
     */
    private Object getValueForInput(Map inputDef, Map complexParams) {
        def paramName = inputDef.name
        def paramType = inputDef.type
        def paramValue = complexParams.get(paramName)
        
        if( paramValue == null ) {
            throw new IllegalArgumentException("Missing required parameter: --${paramName}")
        }
        
        // Type-specific conversion
        switch( paramType ) {
            case 'path':
            case 'file':
                if( paramValue instanceof String || paramValue instanceof GString ) {
                    return parseFileInput(paramValue.toString())
                }
                return paramValue
                
            case 'val':
                // For val inputs, return as-is (could be Map for complex structures)
                return paramValue
                
            case 'env':
                return paramValue?.toString()
                
            default:
                return paramValue
        }
    }
}
