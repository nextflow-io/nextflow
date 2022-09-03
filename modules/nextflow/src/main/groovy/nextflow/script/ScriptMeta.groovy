/*
 * Copyright 2020-2022, Seqera Labs
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

import java.lang.reflect.Method
import java.lang.reflect.Modifier
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.exception.DuplicateModuleFunctionException
import nextflow.exception.MissingModuleComponentException
import nextflow.script.bundle.ResourcesBundle
/**
 * Holds a nextflow script meta-data such as the
 * defines processes and workflows, the included modules
 * the script path, etc.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ScriptMeta {

    static private List<String> INVALID_FUNCTION_NAMES = [
            'methodMissing',
            'propertyMissing'
    ]

    static private Map<BaseScript,ScriptMeta> REGISTRY = new HashMap<>(10)

    static private Set<String> resolvedProcessNames = new HashSet<>(20)

    static ScriptMeta get(BaseScript script) {
        if( !script ) throw new IllegalStateException("Missing current script context")
        return REGISTRY.get(script)
    }

    static Set<String> allProcessNames() {
        def result = new HashSet()
        for( ScriptMeta entry : REGISTRY.values() )
            result.addAll( entry.getProcessNames() )
        // add all resolved names
        result.addAll(resolvedProcessNames)
        return result
    }

    static void addResolvedName(String name) {
        resolvedProcessNames.add(name)
    }

    static Map<String,Path> allScriptNames() {
        def result = new HashMap(REGISTRY.size())
        for( ScriptMeta entry : REGISTRY.values() )
            result.put(entry.scriptName, entry.scriptPath)
        return result
    }

    static ScriptMeta current() {
        get(ExecutionStack.script())
    }

    /** the script {@link Class} object */
    private Class<? extends BaseScript> clazz

    /** The location path from where the script has been loaded */
    private Path scriptPath

    /** The list of function, procs and workflow defined in this script */
    private Map<String,ComponentDef> definitions = new HashMap<>(10)

    /** The module components included in the script */
    private Map<String,ComponentDef> imports = new HashMap<>(10)

    private List<String> dsl1ProcessNames

    /** Whenever it's a module script or the main script */
    private boolean module

    private Map<String,Integer> functionsCount = new HashMap<>()

    Path getScriptPath() { scriptPath }

    Path getModuleDir () { scriptPath?.parent }

    String getScriptName() { clazz.getName() }

    boolean isModule() { module }

    ScriptMeta(BaseScript script) {
        this.clazz = script.class
        for( def entry : definedFunctions0(script) ) {
            addDefinition(entry)
        }
    }

    /** only for testing */
    protected ScriptMeta() {}

    @PackageScope
    void setScriptPath(Path path) {
        scriptPath = path
    }

    @PackageScope
    void setModule(boolean val) {
        this.module = val
    }

    private void incFunctionCount(String name) {
        final count = functionsCount.getOrDefault(name, 0)
        functionsCount.put(name, count+1)
    }

    void validate() {
        // check for duplicate function names
        for( final name : functionsCount.keySet() ) {
            if( functionsCount.get(name)<2 )
                continue
            final msg = "A function with name '$name' is defined more than once in module script: $scriptPath -- Make sure to not define the same function with multiple signatures or arguments with a default value"
            if( NF.isStrictMode() )
                throw new DuplicateModuleFunctionException(msg)
            log.warn(msg)
        }
    }

    void checkComponentName(ComponentDef component, String name) {
        if( component !instanceof ProcessDef && component !instanceof FunctionDef ) {
            return
        }
        if (functionsCount.get(component.name)) {
            final msg = "A function with name '$name' is defined more than once in module script: $scriptPath -- Make sure to not define the same function as process"
            if (NF.isStrictMode())
                throw new DuplicateModuleFunctionException(msg)
            log.warn(msg)
        }
        if (imports.get(component.name)) {
            final msg = "A process with name '$name' is defined more than once in module script: $scriptPath -- Make sure to not define the same function as process"
            if (NF.isStrictMode())
                throw new DuplicateModuleFunctionException(msg)
            log.warn(msg)
        }
    }

    /*
     * This method invocation is made by the NF AST transformer to pass
     * the process names declared in the workflow script. This is only required
     * for DSL1 script.
     *
     * When using DSL2 process names can be discovered during
     * the script execution since, the process declaration is de-coupled by the
     * process invocations.
     */
    @PackageScope
    void setDsl1ProcessNames(List<String> names) {
        this.dsl1ProcessNames = names
    }

    @PackageScope
    List<String> getDsl1ProcessNames() {
        dsl1ProcessNames ?: Collections.<String>emptyList()
    }

    @PackageScope
    static ScriptMeta register(BaseScript script) {
        def meta = new ScriptMeta(script)
        REGISTRY.put(script, meta)
        return meta
    }

    static List<FunctionDef> definedFunctions0(BaseScript script) {
        final allMethods = script.class.getDeclaredMethods()
        final result = new ArrayList<FunctionDef>(allMethods.length)
        for( Method method : allMethods ) {
            if( !Modifier.isPublic(method.getModifiers()) ) continue
            if( Modifier.isStatic(method.getModifiers())) continue
            if( method.name.startsWith('super$')) continue
            if( method.name in INVALID_FUNCTION_NAMES ) continue

            // If method is already into the list, maybe with other signature, it's not necessary to include it again
            if( result.find{it.name == method.name}) continue
            result.add(new FunctionDef(script, method.name))
        }
        return result
    }

    ScriptMeta addDefinition(ComponentDef component) {
        final name = component.name
        if( !module && NF.hasOperator(name) )
            log.warn "${component.type.capitalize()} with name '$name' overrides a built-in operator with the same name"
        checkComponentName(component, name)
        definitions.put(component.name, component)
        if( component instanceof FunctionDef ){
            incFunctionCount(name)
        }
        return this
    }

    ScriptMeta addDefinition(ComponentDef... component) {
        for( def entry : component ) {
            addDefinition(entry)
        }
        return this
    }

    Collection<ComponentDef> getDefinitions() {
        return definitions.values()
    }

    ComponentDef getComponent(String name) {
        definitions.get(name) ?: imports.get(name)
    }

    WorkflowDef getWorkflow(String name) {
        (WorkflowDef)getComponent(name)
    }

    ProcessDef getProcess(String name) {
        (ProcessDef)getComponent(name)
    }

    FunctionDef getFunction(String name) {
        (FunctionDef)getComponent(name)
    }

    Set<String> getAllNames() {
        final result = new HashSet(definitions.size() + imports.size())
        // local definitions
        for( def item : definitions.values() ) {
            if( item.name )
                result.add(item.name)
        }
        // processes from imports
        for( def item: imports.values() ) {
            if( item.name )
                result.add(item.name)
        }
        return result
    }

    Set<String> getProcessNames() {
        if( NF.dsl1 )
            return new HashSet<String>(getDsl1ProcessNames())

        def result = new HashSet(definitions.size() + imports.size())
        // local definitions
        for( def item : definitions.values() ) {
            if( item instanceof ProcessDef )
                result.add(item.name)
        }
        // processes from imports
        for( def item: imports.values() ) {
            if( item instanceof ProcessDef )
                result.add(item.name)
        }
        return result
    }

    Set<String> getLocalProcessNames() {
        if( NF.dsl1 )
            return new HashSet<String>(getDsl1ProcessNames())

        def result = new HashSet(definitions.size() + imports.size())
        // local definitions
        for( def item : definitions.values() ) {
            if( item instanceof ProcessDef )
                result.add(item.name)
        }
        return result
    }

    Set<String> getLocalWorkflowNames() {
        def result = new HashSet(definitions.size())
        for( def item : definitions.values() ) {
            if( item instanceof WorkflowDef && item.name )
                result.add(item.name)
        }
        return result
    }

    void addModule(BaseScript script, String name, String alias) {
       addModule(get(script), name, alias)
    }

    void addModule(ScriptMeta script, String name, String alias) {
        assert script
        if( name ) {
            // include a specific
            def item = script.getComponent(name)
            if( !item )
                throw new MissingModuleComponentException(script, name)
            addModule0(item, alias)
        }
        else for( def item : script.getDefinitions() ) {
            addModule0(item)
        }
    }

    protected void addModule0(ComponentDef component, String alias=null) {
        assert component

        final name = alias ?: component.name
        checkComponentName(component, name)
        if( name != component.name ) {
            imports.put(name, component.cloneWithName(name))
        }
        else {
            imports.put(name, component)
        }
    }

    @Memoized
    ResourcesBundle getModuleBundle() {
        if( !scriptPath )
            throw new IllegalStateException("Module scriptPath has not been defined yet")
        if( scriptPath.getName()!='main.nf' )
            return null
        final bundlePath = scriptPath.resolveSibling('resources')
        return ResourcesBundle.scan(bundlePath)
    }

}
