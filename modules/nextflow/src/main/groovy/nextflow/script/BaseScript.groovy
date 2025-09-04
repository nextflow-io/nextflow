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

import java.lang.reflect.InvocationTargetException
import java.nio.file.Paths

import groovy.util.logging.Slf4j
import nextflow.Channel
import nextflow.NF
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.secret.SecretsLoader

/**
 * Any user defined script will extends this class, it provides the base execution context
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class BaseScript extends Script implements ExecutionContext {

    private Session session

    private ProcessFactory processFactory

    private ScriptMeta meta

    private WorkflowDef entryFlow

    private OutputDef publisher

    @Lazy InputStream stdin = { System.in }()

    BaseScript() {
        meta = ScriptMeta.register(this)
    }

    BaseScript(Binding binding) {
        super(binding)
        meta = ScriptMeta.register(this)
    }

    @Override
    ScriptBinding getBinding() {
        (ScriptBinding)super.getBinding()
    }

    Session getSession() {
        session
    }

    /**
     * Holds the configuration object which will used to execution the user tasks
     */
    @Deprecated
    protected Map getConfig() {
        final msg = "The access of `config` object is deprecated"
        throw new DeprecationException(msg)
    }

    /**
     * Enable disable task 'echo' configuration property
     * @param value
     */
    @Deprecated
    protected void echo(boolean value = true) {
        final msg = "The use of `echo` method has been deprecated"
        throw new DeprecationException(msg)
    }

    private void setup() {
        binding.owner = this
        session = binding.getSession()
        processFactory = session.newProcessFactory(this)

        binding.setVariable( 'baseDir', session.baseDir )
        binding.setVariable( 'projectDir', session.baseDir )
        binding.setVariable( 'workDir', session.workDir )
        binding.setVariable( 'workflow', session.workflowMetadata )
        binding.setVariable( 'nextflow', NextflowMeta.instance )
        binding.setVariable( 'launchDir', Paths.get('./').toRealPath() )
        binding.setVariable( 'moduleDir', meta.moduleDir )
        binding.setVariable( 'secrets', SecretsLoader.secretContext() )
    }

    protected void params(Closure body) {
        if( entryFlow )
            throw new IllegalStateException("Workflow params definition must be defined before the entry workflow")
        if( ExecutionStack.withinWorkflow() )
            throw new IllegalStateException("Workflow params definition is not allowed within a workflow")

        final dsl = new ParamsDsl()
        final cl = (Closure)body.clone()
        cl.setDelegate(dsl)
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.call()

        dsl.apply(session)
    }

    protected process( String name, Closure<BodyDef> body ) {
        final process = new ProcessDef(this,body,name)
        meta.addDefinition(process)
    }

    /**
     * Workflow main entry point
     *
     * @param body The implementation body of the workflow
     * @return The result of workflow execution
     */
    protected workflow(Closure<BodyDef> workflowBody) {
        // launch the execution
        final workflow = new WorkflowDef(this, workflowBody)
        // capture the main (unnamed) workflow definition
        this.entryFlow = workflow
        // add it to the list of workflow definitions
        meta.addDefinition(workflow)
    }

    protected workflow(String name, Closure<BodyDef> workflowDef) {
        final workflow = new WorkflowDef(this,workflowDef,name)
        meta.addDefinition(workflow)
    }

    protected output(Closure closure) {
        if( !NF.outputDefinitionEnabled )
            throw new IllegalStateException("Workflow output definition requires the `nextflow.preview.output` feature flag")
        if( !entryFlow )
            throw new IllegalStateException("Workflow output definition must be defined after the entry workflow")
        if( ExecutionStack.withinWorkflow() )
            throw new IllegalStateException("Workflow output definition is not allowed within a workflow")

        publisher = new OutputDef(closure)
    }

    protected IncludeDef include( IncludeDef include ) {
        if(ExecutionStack.withinWorkflow())
            throw new IllegalStateException("Include statement is not allowed within a workflow definition")
        include .setSession(session)
    }

    /**
     * Invokes custom methods in the task execution context
     *
     * @see nextflow.processor.TaskContext#invokeMethod(java.lang.String, java.lang.Object)
     * @see WorkflowBinding#invokeMethod(java.lang.String, java.lang.Object)
     *
     * @param name the name of the method to call
     * @param args the arguments to use for the method call
     * @return The result of the custom method execution
     */
    @Override
    Object invokeMethod(String name, Object args) {
        binding.invokeMethod(name, args)
    }

    private run0() {
        final result = runScript()
        if( meta.isModule() ) {
            return result
        }

        // if an `entryName` was specified via the command line, resolve it to an entryFlow
        if( binding.entryName ) {
            // Check for process entry syntax: 'process:NAME'
            if( binding.entryName.startsWith('process:') ) {
                final processName = binding.entryName.substring(8) // Remove 'process:' prefix
                final processDef = meta.getProcess(processName)
                if( !processDef ) {
                    def msg = "Unknown process entry name: ${processName}"
                    final allProcessNames = meta.getProcessNames()
                    final guess = allProcessNames.closest(processName)
                    if( guess )
                        msg += " -- Did you mean?\n" + guess.collect { "  $it"}.join('\n')
                    throw new IllegalArgumentException(msg)
                }
                // Create a workflow to execute the specified process with parameter mapping
                entryFlow = createProcessEntryWorkflow(processDef)
            }
            // Traditional workflow entry
            else if( !(entryFlow=meta.getWorkflow(binding.entryName) ) ) {
                def msg = "Unknown workflow entry name: ${binding.entryName}"
                final allNames = meta.getWorkflowNames()
                final guess = allNames.closest(binding.entryName)
                if( guess )
                    msg += " -- Did you mean?\n" + guess.collect { "  $it"}.join('\n')
                throw new IllegalArgumentException(msg)
            }
        }

        if( !entryFlow ) {
            if( meta.getLocalWorkflowNames() )
                throw new AbortOperationException("No entry workflow specified")
            // Check if we have a single standalone process that can be executed automatically
            if( meta.hasSingleExecutableProcess() ) {
                // Create a workflow to execute the single process
                entryFlow = createSingleProcessWorkflow()
            }
            // Check if we have multiple processes that require -entry specification
            else if( meta.hasMultipleExecutableProcesses() ) {
                def processNames = meta.getLocalProcessNames()
                throw new AbortOperationException("Multiple processes found (${processNames.join(', ')}). Use -entry process:NAME to specify which process to execute.")
            } else {
                return result
            }
        }

        // invoke the entry workflow
        session.notifyBeforeWorkflowExecution()
        final ret = entryFlow.invoke_a(BaseScriptConsts.EMPTY_ARGS)
        if( publisher )
            publisher.apply(session)
        session.notifyAfterWorkflowExecution()
        return ret
    }

    Object run() {
        setup()
        ExecutionStack.push(this)
        try {
            run0()
        }
        catch( InvocationTargetException e ) {
            // provide the exception cause which is more informative than InvocationTargetException
            Throwable target = e
            do target = target.cause
            while ( target instanceof InvocationTargetException )
            throw target
        }
        finally {
            ExecutionStack.pop()
        }
    }

    protected abstract Object runScript()

    //========================================
    // PROCESS ENTRY EXECUTION FEATURE
    //========================================
    /*
     * The methods below implement the process entry execution feature, which allows
     * users to execute Nextflow processes directly without writing explicit workflows.
     * 
     * Key functionality:
     * - Single process scripts run automatically: `nextflow run script.nf --param value`
     * - Multi-process scripts use entry selection: `nextflow run script.nf -entry process:name --param value`
     * - Command-line parameters are automatically mapped to process input channels
     * - Supports all standard Nextflow input types: val, path, env, tuple, each
     * 
     * Implementation approach:
     * 1. Parse the process body to extract input parameter definitions
     * 2. Map command-line arguments to the extracted parameter specifications
     * 3. Create appropriate Nextflow channels for each mapped parameter
     * 4. Generate synthetic workflows that execute the process with mapped inputs
     * 
     * This feature bridges the gap between command-line tools and workflow orchestration,
     * making Nextflow processes more accessible for direct execution scenarios.
     */

    /**
     * Creates a workflow to execute a single standalone process automatically.
     * This allows single-process scripts to run without requiring the -entry option.
     * 
     * The method validates that exactly one process exists and creates a synthetic
     * workflow that maps command-line parameters to the process inputs.
     *
     * @return WorkflowDef that executes the single process with parameter mapping
     * @throws IllegalStateException if there is not exactly one process
     */
    protected WorkflowDef createSingleProcessWorkflow() {
        // Validate preconditions - must have exactly one process
        def processNames = meta.getLocalProcessNames()
        if( processNames.size() != 1 ) {
            throw new IllegalStateException("Expected exactly one process, found: ${processNames.size()}")
        }
        
        final processName = processNames.first()
        final processDef = meta.getProcess(processName)
        
        // Create the workflow execution logic for the single process
        def workflowLogic = { ->
            log.debug "Executing auto-generated single process workflow for: ${processName}"
            
            // Map command-line parameters to process input channels
            def inputChannels = createProcessInputChannelsWithMapping(processDef)
            
            log.debug "Mapped ${inputChannels.length} input channels for single process: ${processName}"
            
            // Execute the single process with the mapped input channels
            this.invokeMethod(processName, inputChannels)
        }
        
        // Create workflow metadata for debugging and introspection
        def sourceCode = "    // Auto-generated single process workflow\n    ${processName}(...)"
        
        // Wrap the logic in a BodyDef closure as expected by WorkflowDef constructor
        def workflowBody = { ->
            return new BodyDef(workflowLogic, sourceCode, 'workflow')
        }
        
        return new WorkflowDef(this, workflowBody)
    }


    /**
     * Creates a workflow to execute a specific process with parameter mapping.
     * This enables process execution via the -entry process:NAME syntax.
     * 
     * The workflow automatically maps command-line parameters to process inputs
     * by analyzing the process definition and creating appropriate channels.
     *
     * @param processDef The ProcessDef object for the target process
     * @return WorkflowDef that executes the process with parameter mapping
     */
    protected WorkflowDef createProcessEntryWorkflow(ProcessDef processDef) {
        final processName = processDef.name
        
        // Create the workflow execution logic that handles parameter mapping
        def workflowLogic = { ->
            log.debug "Executing process entry workflow for: ${processName}"
            
            // Map command-line parameters to process input channels
            def inputChannels = createProcessInputChannelsWithMapping(processDef)
            
            log.debug "Mapped ${inputChannels.length} input channels for process: ${processName}"
            
            // Execute the process with the mapped input channels
            this.invokeMethod(processName, inputChannels)
        }
        
        // Create workflow metadata for debugging and introspection
        def sourceCode = "    // Auto-generated process entry workflow\n    ${processName}(...)"
        
        // Wrap the logic in a BodyDef closure as expected by WorkflowDef constructor
        def workflowBody = { ->
            return new BodyDef(workflowLogic, sourceCode, 'workflow')
        }
        
        return new WorkflowDef(this, workflowBody)
    }

    /**
     * Creates input channels for a process by mapping command-line parameters to process inputs.
     * 
     * This is the main parameter mapping method that:
     * 1. Extracts input definitions from the process body by parsing Nextflow DSL
     * 2. Maps command-line parameters to the extracted input specifications  
     * 3. Creates appropriate Nextflow channels for each input based on its type
     * 4. Handles missing parameters with appropriate error messages
     *
     * @param processDef The ProcessDef object containing the process definition
     * @return Object[] Array of channels corresponding to process inputs, empty if no inputs or error
     */
    protected Object[] createProcessInputChannelsWithMapping(ProcessDef processDef) {
        final processName = processDef.name
        log.debug "Starting parameter mapping for process: ${processName}"
        
        try {
            // Step 1: Extract input definitions from process body
            List<Map> inputSpecs = extractProcessInputDefinitions(processDef)
            
            if( inputSpecs.isEmpty() ) {
                log.debug "No input parameters found for process: ${processName}"
                return new Object[0]
            }
            
            log.debug "Found ${inputSpecs.size()} input parameters for process ${processName}: ${inputSpecs.collect { it.name }}"
            
            // Step 2: Map command-line parameters to channels
            Object[] inputChannels = mapParametersToChannels(inputSpecs)
            
            log.debug "Successfully mapped ${inputChannels.length} input channels for process: ${processName}"
            return inputChannels
            
        } catch (Exception e) {
            log.warn "Failed to create input channels for process ${processName}: ${e.message}"
            return new Object[0]
        }
    }
    
    /**
     * Extracts input parameter definitions from a process body by parsing the Nextflow DSL.
     * 
     * This method uses a specialized delegate to intercept internal Nextflow method calls
     * (_in_val, _in_path, etc.) that represent input declarations in the compiled process body.
     * 
     * @param processDef The ProcessDef to analyze
     * @return List of Maps containing input specifications [type: String, name: String]
     */
    private List<Map> extractProcessInputDefinitions(ProcessDef processDef) {
        // Access the raw process body closure - this contains the compiled Nextflow DSL
        def rawBody = processDef.rawBody
        if( !rawBody ) {
            log.debug "No process body found for: ${processDef.name}"
            return []
        }
        
        log.debug "Analyzing process body for input definitions: ${processDef.name}"
        
        // Clone the body to avoid side effects and set up input extraction
        def bodyClone = rawBody.clone()
        def extractionDelegate = new ProcessInputExtractionDelegate()
        bodyClone.setDelegate(extractionDelegate)
        bodyClone.setResolveStrategy(Closure.DELEGATE_FIRST)
        
        try {
            // Execute the cloned body - this will trigger our delegate methods
            bodyClone.call()
            
            // Return the collected input definitions
            def extractedInputs = extractionDelegate.getExtractedInputs()
            log.debug "Extracted ${extractedInputs.size()} input definitions from process ${processDef.name}"
            return extractedInputs
            
        } catch (Exception e) {
            log.debug "Failed to extract input definitions from process ${processDef.name}: ${e.message}"
            return []
        }
    }
    
    /**
     * Maps input parameter specifications to Nextflow channels using command-line parameter values.
     * 
     * For each input specification, this method:
     * 1. Looks up the parameter value in session.params (populated from command-line args)
     * 2. Creates an appropriate Nextflow channel based on the input type (val, path, etc.)
     * 3. Handles missing required parameters with descriptive error messages
     *
     * @param inputSpecs List of input specifications from process definition
     * @return Object[] Array of Nextflow channels for process execution
     */
    private Object[] mapParametersToChannels(List<Map> inputSpecs) {
        Object[] channels = new Object[inputSpecs.size()]
        
        log.debug "Available command-line parameters: ${session.params.keySet()}"
        
        // Use traditional for loop for better performance and clearer index handling
        for( int i = 0; i < inputSpecs.size(); i++ ) {
            Map inputSpec = inputSpecs[i]
            String paramName = inputSpec.name
            String inputType = inputSpec.type
            
            log.debug "Processing parameter '${paramName}' of type '${inputType}'"
            
            // Look up parameter value from command-line arguments
            def paramValue = session.params.get(paramName)
            
            if( paramValue != null ) {
                log.debug "Found parameter value: ${paramName} = ${paramValue}"
                channels[i] = createChannelForInputType(inputType, paramValue)
            } else {
                log.debug "Parameter '${paramName}' not provided via command-line"
                channels[i] = createDefaultChannelForInputType(inputType, paramName)
            }
        }
        
        return channels
    }

    /**
     * Creates a Nextflow channel for a process input parameter based on its type and value.
     * 
     * This method handles the conversion from command-line parameter values to the appropriate
     * Nextflow channel types required by different input parameter declarations.
     *
     * @param inputType The type of input parameter (val, path, env, tuple, each)
     * @param paramValue The parameter value from command-line arguments
     * @return Nextflow channel containing the parameter value(s)
     */
    protected Object createChannelForInputType(String inputType, def paramValue) {
        switch( inputType ) {
            case 'val':
                // Value parameters: direct channel creation with the parameter value
                return Channel.of(paramValue)
                
            case 'path':
                // Path parameters: convert strings to Path objects and validate existence
                def path = (paramValue instanceof String) ? 
                    java.nio.file.Paths.get(paramValue) : paramValue
                    
                // Warn if file doesn't exist (non-fatal for flexibility)
                if( path instanceof java.nio.file.Path && !java.nio.file.Files.exists(path) ) {
                    log.warn "Path parameter references non-existent file: ${path}"
                }
                return Channel.of(path)
                
            case 'env':
                // Environment parameters: pass through as value channel
                return Channel.of(paramValue)
                
            case 'tuple':
                // Tuple parameters: handle composite values
                if( paramValue instanceof Collection ) {
                    return Channel.of(paramValue as List)
                } else {
                    // Wrap single values in a list for tuple consistency
                    return Channel.of([paramValue])
                }
                
            case 'each':
                // Each parameters: create iterable channels for collection processing
                if( paramValue instanceof Collection ) {
                    return Channel.fromIterable(paramValue)
                } else if( paramValue instanceof String && paramValue.contains(',') ) {
                    // Handle comma-separated values as collections
                    def items = paramValue.split(',').collect { it.trim() }
                    return Channel.fromIterable(items)
                } else {
                    return Channel.of(paramValue)
                }
                
            default:
                // Unknown input types: default to value channel with warning
                log.debug "Unknown input parameter type '${inputType}', using value channel"
                return Channel.of(paramValue)
        }
    }

    /**
     * Creates a default channel or throws an error when a required parameter is missing.
     * 
     * This method handles missing command-line parameters by either providing sensible
     * defaults (for optional parameters like env) or throwing descriptive error messages
     * for required parameters.
     *
     * @param inputType The type of input parameter that is missing
     * @param paramName The name of the missing parameter for error messages
     * @return Default channel for optional parameters
     * @throws IllegalArgumentException for required parameters that are missing
     */
    protected Object createDefaultChannelForInputType(String inputType, String paramName) {
        switch( inputType ) {
            case 'val':
                // Value parameters are typically required
                throw new IllegalArgumentException("Missing required value parameter: --${paramName}")
                
            case 'path':
                // Path parameters are typically required 
                throw new IllegalArgumentException("Missing required path parameter: --${paramName}")
                
            case 'env':
                // Environment parameters may have defaults or be optional
                // Return empty channel to allow process to handle gracefully
                return Channel.empty()
                
            case 'tuple':
                // Tuple parameters are typically required
                throw new IllegalArgumentException("Missing required tuple parameter: --${paramName}")
                
            case 'each':
                // Each parameters are typically required for iteration
                throw new IllegalArgumentException("Missing required each parameter: --${paramName}")
                
            default:
                // Unknown parameter types: assume required and provide generic error
                throw new IllegalArgumentException("Missing required parameter: --${paramName}")
        }
    }

    /**
     * Specialized delegate class for extracting input parameter definitions from Nextflow process bodies.
     * 
     * This class intercepts method calls when a cloned process body is executed, allowing us to
     * capture input parameter declarations without actually executing the full process logic.
     * 
     * The Nextflow compiler converts input declarations like "val name" into internal method calls
     * like "_in_val(TokenVar(name))", which this delegate captures and converts back to 
     * structured parameter specifications.
     */
    protected class ProcessInputExtractionDelegate {
        
        /** List of extracted input specifications in format [type: String, name: String] */
        private final List<Map> extractedInputs = []
        
        /**
         * Returns the collected input parameter definitions.
         * @return List of input specifications
         */
        List<Map> getExtractedInputs() {
            return new ArrayList<>(extractedInputs)
        }
        
        /**
         * Handles traditional input {} block declarations (less common in compiled code).
         * @param inputClosure Closure containing input declarations
         */
        def input(Closure inputClosure) {
            log.debug "Processing traditional input block"
            def inputDelegate = new TraditionalInputParsingDelegate()
            inputClosure.setDelegate(inputDelegate)
            inputClosure.setResolveStrategy(Closure.DELEGATE_FIRST)
            inputClosure.call()
            extractedInputs.addAll(inputDelegate.getInputs())
        }
        
        /*
         * Methods below handle internal Nextflow DSL method calls that represent
         * compiled input declarations. These are generated by the Nextflow compiler
         * from user DSL like "val name" -> "_in_val(TokenVar(name))".
         */
        
        /** Handles value input parameters: val paramName */
        def _in_val(Object tokenVar) {
            addInputParameter('val', tokenVar)
        }
        
        /** Handles file input parameters: file paramName (legacy syntax) */
        def _in_file(Object tokenVar) {
            addInputParameter('path', tokenVar)
        }
        
        /** Handles path input parameters: path paramName */
        def _in_path(Object tokenVar) {
            addInputParameter('path', tokenVar)
        }
        
        /** Handles environment input parameters: env paramName */
        def _in_env(Object tokenVar) {
            addInputParameter('env', tokenVar)
        }
        
        /** Handles each input parameters: each paramName */
        def _in_each(Object tokenVar) {
            addInputParameter('each', tokenVar)
        }
        
        /** Handles tuple input parameters: tuple paramName1, paramName2, ... */
        def _in_tuple(Object... tokenVars) {
            log.debug "Processing tuple input with ${tokenVars.length} elements"
            
            // Process each element of the tuple as a separate parameter
            for( int i = 0; i < tokenVars.length; i++ ) {
                addInputParameter('tuple', tokenVars[i])
            }
        }
        
        /**
         * Common helper method to extract parameter name from TokenVar and add to collection.
         * @param inputType The type of input parameter (val, path, etc.)
         * @param tokenVar The TokenVar object containing the parameter name
         */
        private void addInputParameter(String inputType, Object tokenVar) {
            if( tokenVar?.hasProperty('name') ) {
                String paramName = tokenVar.name
                log.debug "Extracted ${inputType} parameter: ${paramName}"
                extractedInputs.add([type: inputType, name: paramName])
            } else {
                log.debug "Skipped ${inputType} parameter with invalid TokenVar: ${tokenVar}"
            }
        }
        
        /**
         * Catches all other method calls and ignores them.
         * This allows the delegate to skip process elements like output, script, etc.
         * @param methodName Name of the method being called
         * @param args Arguments passed to the method
         */
        def methodMissing(String methodName, Object args) {
            log.debug "Ignoring process element: ${methodName}"
            // Silently ignore other process declarations (output, script, etc.)
        }
    }

    /**
     * Delegate class for parsing traditional input {} block declarations.
     * 
     * This handles the less common case where users write explicit input blocks
     * in their process definitions. Most modern Nextflow code uses the direct
     * syntax (e.g., "val name" instead of "input { val name }").
     */
    protected class TraditionalInputParsingDelegate {
        
        /** List of parsed input specifications */
        private final List<Map> inputs = []
        
        /**
         * Returns the collected input specifications.
         * @return List of input specifications
         */
        List<Map> getInputs() {
            return new ArrayList<>(inputs)
        }
        
        /** Handles value input declarations: val name */
        def val(String name) {
            inputs.add([type: 'val', name: name])
        }
        
        /** Handles path input declarations: path name */
        def path(String name) {
            inputs.add([type: 'path', name: name])
        }
        
        /** Handles file input declarations: file name (legacy) */
        def file(String name) {
            inputs.add([type: 'path', name: name])
        }
        
        /** Handles environment input declarations: env name */
        def env(String name) {
            inputs.add([type: 'env', name: name])
        }
        
        /** Handles tuple input declarations: tuple name */
        def tuple(String name) {
            inputs.add([type: 'tuple', name: name])
        }
        
        /** Handles each input declarations: each name */
        def each(String name) {
            inputs.add([type: 'each', name: name])
        }
        
        /**
         * Generic handler for any other input types that might be added in the future.
         * @param inputType The name of the input method called
         * @param args Arguments passed to the method
         */
        def methodMissing(String inputType, Object args) {
            log.debug "Processing generic input type: ${inputType} with args: ${args}"
            
            if( args instanceof Object[] && args.length > 0 ) {
                String name = args[0]?.toString()
                if( name ) {
                    log.debug "Added generic input: [type: ${inputType}, name: ${name}]"
                    inputs.add([type: inputType, name: name])
                } else {
                    log.debug "Skipped input with invalid name: ${inputType}"
                }
            }
        }
    }

    @Override
    void print(Object object) {
        if( session?.quiet )
            return

        if( session?.ansiLog )
            log.info(object?.toString())
        else
            super.print(object)
    }

    @Override
    void println() {
        if( session?.quiet )
            return

        if( session?.ansiLog )
            log.info("")
        else
            super.println()
    }

    @Override
    void println(Object object) {
        if( session?.quiet )
            return

        if( session?.ansiLog )
            log.info(object?.toString())
        else
            super.println(object)
    }

    @Override
    void printf(String msg, Object arg) {
        if( session?.quiet )
            return

        if( session?.ansiLog )
            log.info(String.printf(msg, arg))
        else
            super.printf(msg, arg)
    }

    @Override
    void printf(String msg, Object[] args) {
        if( session?.quiet )
            return

        if( session?.ansiLog )
            log.info(String.printf(msg, args))
        else
            super.printf(msg, args)
    }

}
