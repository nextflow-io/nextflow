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

import java.lang.reflect.Method
import java.lang.reflect.Modifier
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.DuplicateModuleIncludeException
import nextflow.exception.MissingModuleComponentException
import nextflow.extension.OperatorEx
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

    static ScriptMeta get(BaseScript script) {
        if( !script ) throw new IllegalStateException("Missing current script context")
        REGISTRY.get(script)
    }

    static Map<String,Path> allScriptNames() {
        def result = new HashMap(REGISTRY.size())
        for( Map.Entry<BaseScript,ScriptMeta> entry : REGISTRY )
            result.put(entry.key.getClass().getName(), entry.value.scriptPath)
        return result
    }

    static ScriptMeta current() {
        get(ExecutionStack.script())
    }

    /** the script {@link Class} object */
    private Class<? extends BaseScript> clazz

    /** The location path from there the script has been loaded */
    private Path scriptPath

    /** The list of function, procs and workflow defined in this script */
    private Map<String,ComponentDef> definitions = new HashMap<>(10)

    /** The module components included in the script */
    private Map<String,ComponentDef> imports = new HashMap<>(10)

    /** Whenever it's a module script or the main script */
    private boolean module

    Path getScriptPath() { scriptPath }

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

    @PackageScope
    static ScriptMeta register(BaseScript script) {
        def meta = new ScriptMeta(script)
        REGISTRY.put(script, meta)
        return meta
    }

    static List<FunctionDef> definedFunctions0(BaseScript script) {
        def allMethods = script.class.getDeclaredMethods()
        def result = new ArrayList(allMethods.length)
        for( Method method : allMethods ) {
            if( !Modifier.isPublic(method.getModifiers()) ) continue
            if( Modifier.isStatic(method.getModifiers())) continue
            if( method.name.startsWith('super$')) continue
            if( method.name in INVALID_FUNCTION_NAMES ) continue

            result.add(new FunctionDef(script, method))
        }
        return result
    }

    ScriptMeta addDefinition(ComponentDef component) {
        final name = component.name
        if( !module && name in OperatorEx.OPERATOR_NAMES )
            log.warn "${component.type.capitalize()} with name '$name' overrides a built-in operator with the same name"
        definitions.put(component.name, component)
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

    Set<String> getProcessNames() {
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
        final existing = getComponent(name)
        if (existing) {
            def msg = "A ${existing.type} with name '$name' is already defined in the current context"
            throw new DuplicateModuleIncludeException(msg)
        }

        if( name != component.name ) {
            imports.put(name, component.withName(name))
        }
        else {
            imports.put(name, component)
        }
    }

}
