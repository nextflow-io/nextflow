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
     * - Command-line parameters are passed directly to processes
     * - Supports all standard Nextflow input types: val, path, env, tuple, each
     * 
     * Implementation approach:
     * 1. Parse the process body to extract input parameter names
     * 2. Look up corresponding values from command-line arguments (session.params)
     * 3. Pass parameter values directly to the process
     * 4. Let Nextflow handle channel creation automatically
     * 
     * This simplified approach relies on Nextflow's built-in parameter handling rather
     * than manually creating channels, making the implementation much cleaner.
     */

    /**
     * Creates a workflow to execute a single standalone process automatically.
     * This allows single-process scripts to run without requiring the -entry option.
     * 
     * The method validates that exactly one process exists and creates a synthetic
     * workflow that passes command-line parameters directly to the process.
     *
     * @return WorkflowDef that executes the single process with direct parameter passing
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
            
            // Get input parameter names and pass corresponding values from session.params
            def inputArgs = getProcessInputArguments(processDef)
            
            log.debug "Passing ${inputArgs.size()} arguments to single process: ${processName}"
            
            // Execute the single process with the parameter values
            this.invokeMethod(processName, inputArgs as Object[])
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
     * Creates a workflow to execute a specific process with direct parameter passing.
     * This enables process execution via the -entry process:NAME syntax.
     * 
     * The workflow passes command-line parameters directly to the process,
     * letting Nextflow handle the channel creation automatically.
     *
     * @param processDef The ProcessDef object for the target process
     * @return WorkflowDef that executes the process with direct parameter passing
     */
    protected WorkflowDef createProcessEntryWorkflow(ProcessDef processDef) {
        final processName = processDef.name
        
        // Create the workflow execution logic that handles parameter passing
        def workflowLogic = { ->
            log.debug "Executing process entry workflow for: ${processName}"
            
            // Get input parameter names and pass corresponding values from session.params
            def inputArgs = getProcessInputArguments(processDef)
            
            log.debug "Passing ${inputArgs.size()} arguments to process: ${processName}"
            
            // Execute the process with the parameter values
            this.invokeMethod(processName, inputArgs as Object[])
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
     * Gets the input arguments for a process by extracting parameter names and looking up
     * their values in session.params (populated from command-line arguments).
     * 
     * This simplified approach lets Nextflow handle channel creation automatically
     * when the process is invoked, rather than manually creating channels.
     * For path parameters, it uses the file() helper to convert strings to Path objects.
     *
     * @param processDef The ProcessDef object containing the process definition
     * @return List of parameter values to pass to the process
     */
    protected List getProcessInputArguments(ProcessDef processDef) {
        final processName = processDef.name
        log.debug "Getting input arguments for process: ${processName}"
        
        try {
            // Extract input parameter names and types from the process body
            List<Map> inputSpecs = extractProcessInputSpecs(processDef)
            
            if( inputSpecs.isEmpty() ) {
                log.debug "No input parameters found for process: ${processName}"
                return []
            }
            
            log.debug "Found ${inputSpecs.size()} input parameters for process ${processName}: ${inputSpecs.collect { it.name }}"
            
            // Get corresponding values from session.params and apply type conversion if needed
            List inputArgs = []
            for( Map inputSpec : inputSpecs ) {
                String paramName = inputSpec.name
                String inputType = inputSpec.type
                def paramValue = session.params.get(paramName)
                
                if( paramValue != null ) {
                    log.debug "Found parameter value: ${paramName} = ${paramValue}"
                    
                    // Convert path parameters using file() helper
                    if( inputType == 'path' && paramValue instanceof String ) {
                        paramValue = nextflow.Nextflow.file(paramValue)
                        log.debug "Converted path parameter ${paramName} to: ${paramValue}"
                    }
                    
                    inputArgs.add(paramValue)
                } else {
                    throw new IllegalArgumentException("Missing required parameter: --${paramName}")
                }
            }
            
            log.debug "Successfully prepared ${inputArgs.size()} arguments for process: ${processName}"
            return inputArgs
            
        } catch (Exception e) {
            log.warn "Failed to get input arguments for process ${processName}: ${e.message}"
            throw e
        }
    }
    
    /**
     * Extracts input parameter specifications (name and type) from a process body by parsing the Nextflow DSL.
     * 
     * This method uses a specialized delegate to intercept internal Nextflow method calls
     * (_in_val, _in_path, etc.) and extracts both parameter names and types.
     * 
     * @param processDef The ProcessDef to analyze
     * @return List of Maps containing [name: String, type: String] in declaration order
     */
    private List<Map> extractProcessInputSpecs(ProcessDef processDef) {
        // Access the raw process body closure - this contains the compiled Nextflow DSL
        def rawBody = processDef.rawBody
        if( !rawBody ) {
            log.debug "No process body found for: ${processDef.name}"
            return []
        }
        
        log.debug "Analyzing process body for input parameter specs: ${processDef.name}"
        
        // Clone the body to avoid side effects and set up input extraction
        def bodyClone = rawBody.clone()
        def extractionDelegate = new ProcessInputSpecExtractor()
        bodyClone.setDelegate(extractionDelegate)
        bodyClone.setResolveStrategy(Closure.DELEGATE_FIRST)
        
        try {
            // Execute the cloned body - this will trigger our delegate methods
            bodyClone.call()
            
            // Return the collected parameter specifications
            def inputSpecs = extractionDelegate.getInputSpecs()
            log.debug "Extracted ${inputSpecs.size()} parameter specs from process ${processDef.name}"
            return inputSpecs
            
        } catch (Exception e) {
            log.debug "Failed to extract parameter specs from process ${processDef.name}: ${e.message}"
            return []
        }
    }
    

    /**
     * Delegate class for extracting parameter names and types from Nextflow process bodies.
     * 
     * This class intercepts method calls when a cloned process body is executed, allowing us to
     * capture input parameter specifications including both names and types.
     * 
     * The Nextflow compiler converts input declarations like "val name" into internal method calls
     * like "_in_val(TokenVar(name))", which this delegate captures.
     */
    protected class ProcessInputSpecExtractor {
        
        /** List of extracted parameter specifications in declaration order */
        private final List<Map> inputSpecs = []
        
        /**
         * Returns the collected parameter specifications.
         * @return List of Maps with [name: String, type: String]
         */
        List<Map> getInputSpecs() {
            return new ArrayList<>(inputSpecs)
        }
        
        /**
         * Handles traditional input {} block declarations (less common in compiled code).
         * @param inputClosure Closure containing input declarations
         */
        def input(Closure inputClosure) {
            log.debug "Processing traditional input block"
            def inputDelegate = new TraditionalInputSpecExtractor()
            inputClosure.setDelegate(inputDelegate)
            inputClosure.setResolveStrategy(Closure.DELEGATE_FIRST)
            inputClosure.call()
            inputSpecs.addAll(inputDelegate.getInputSpecs())
        }
        
        /*
         * Methods below handle internal Nextflow DSL method calls that represent
         * compiled input declarations. These are generated by the Nextflow compiler
         * from user DSL like "val name" -> "_in_val(TokenVar(name))".
         */
        
        /** Handles value input parameters: val paramName */
        def _in_val(Object tokenVar) {
            addInputSpec('val', tokenVar)
        }
        
        /** Handles file input parameters: file paramName (legacy syntax) */
        def _in_file(Object tokenVar) {
            addInputSpec('path', tokenVar)
        }
        
        /** Handles path input parameters: path paramName */
        def _in_path(Object tokenVar) {
            addInputSpec('path', tokenVar)
        }
        
        /** Handles environment input parameters: env paramName */
        def _in_env(Object tokenVar) {
            addInputSpec('env', tokenVar)
        }
        
        /** Handles each input parameters: each paramName */
        def _in_each(Object tokenVar) {
            addInputSpec('each', tokenVar)
        }
        
        /** Handles tuple input parameters: tuple paramName1, paramName2, ... */
        def _in_tuple(Object... tokenVars) {
            log.debug "Processing tuple input with ${tokenVars.length} elements"
            
            // Process each element of the tuple
            for( int i = 0; i < tokenVars.length; i++ ) {
                addInputSpec('tuple', tokenVars[i])
            }
        }
        
        /**
         * Common helper method to extract parameter name from TokenVar and add specification.
         * @param inputType The type of input parameter (val, path, etc.)
         * @param tokenVar The TokenVar object containing the parameter name
         */
        private void addInputSpec(String inputType, Object tokenVar) {
            if( tokenVar?.hasProperty('name') ) {
                String paramName = tokenVar.name
                log.debug "Extracted ${inputType} parameter: ${paramName}"
                inputSpecs.add([name: paramName, type: inputType])
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
     * in their process definitions. Extracts both parameter names and types.
     */
    protected class TraditionalInputSpecExtractor {
        
        /** List of parsed parameter specifications */
        private final List<Map> inputSpecs = []
        
        /**
         * Returns the collected parameter specifications.
         * @return List of Maps with [name: String, type: String]
         */
        List<Map> getInputSpecs() {
            return new ArrayList<>(inputSpecs)
        }
        
        /** Handles value input declarations: val name */
        def val(String name) {
            inputSpecs.add([name: name, type: 'val'])
        }
        
        /** Handles path input declarations: path name */
        def path(String name) {
            inputSpecs.add([name: name, type: 'path'])
        }
        
        /** Handles file input declarations: file name (legacy) */
        def file(String name) {
            inputSpecs.add([name: name, type: 'path'])
        }
        
        /** Handles environment input declarations: env name */
        def env(String name) {
            inputSpecs.add([name: name, type: 'env'])
        }
        
        /** Handles tuple input declarations: tuple name */
        def tuple(String name) {
            inputSpecs.add([name: name, type: 'tuple'])
        }
        
        /** Handles each input declarations: each name */
        def each(String name) {
            inputSpecs.add([name: name, type: 'each'])
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
                    log.debug "Added generic parameter spec: [name: ${name}, type: ${inputType}]"
                    inputSpecs.add([name: name, type: inputType])
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
