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

package nextflow.config

import java.nio.file.Path

import ch.grengine.Grengine
import com.google.common.hash.Hashing
import groovy.transform.PackageScope
import nextflow.ast.NextflowXform
import nextflow.file.FileHelper
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.control.customizers.ImportCustomizer
import org.codehaus.groovy.runtime.InvokerHelper

/*
 * Copyright 2003-2013 the original author or authors.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */



/**
 * A ConfigSlurper that allows to include a config file into another. For example:
 *
 * <pre>
 *     process {
 *         foo = 1
 *         that = 2
 *
 *         includeConfig( 'path/to/another/config/file' )
 *
 *     }
 *
 * </pre>
 *
 * See http://naleid.com/blog/2009/07/30/modularizing-groovy-config-files-with-a-dash-of-meta-programming
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 *
 */

/**
 * <p>
 * ConfigSlurper is a utility class for reading configuration files defined in the form of Groovy
 * scripts. Configuration settings can be defined using dot notation or scoped using closures
 *
 * <pre><code>
 *   grails.webflow.stateless = true
 *    smtp {
 *        mail.host = 'smtp.myisp.com'
 *        mail.auth.user = 'server'
 *    }
 *    resources.URL = "http://localhost:80/resources"
 * </pre></code>
 *
 * <p>Settings can either be bound into nested maps or onto a specified JavaBean instance. In the case
 * of the latter an error will be thrown if a property cannot be bound.
 *
 * @author Graeme Rocher
 * @author Andres Almiray
 * @since 1.5
 */
class ConfigParser {
    private static final ENVIRONMENTS_METHOD = 'environments'

    private Map bindingVars = [:]
    private Map paramVars = [:]

    private final Map<String, List<String>> conditionValues = [:]
    private final Stack<Map<String, ConfigObject>> conditionalBlocks = new Stack<Map<String,ConfigObject>>()
    private final Set<String> conditionalNames = new HashSet<>()
    private final Set<String> profileNames = new HashSet<>()

    private boolean ignoreIncludes

    private boolean renderClosureAsString

    private Grengine grengine

    ConfigParser() {
        this('')
    }

    /**
     * Constructs a new IncludeConfigSlurper instance using the given environment
     * @param env The Environment to use
     */
    ConfigParser(String env) {
        conditionValues[ENVIRONMENTS_METHOD] = [env]
    }

    ConfigParser registerConditionalBlock(String blockName, String blockValue) {
        if (blockName) {
            if (!blockValue) {
                conditionValues.remove(blockName)
            }
            else {
                conditionValues[blockName] = [blockValue]
            }
        }
        return this
    }

    ConfigParser registerConditionalBlock(String blockName, List<String> blockValues) {
        if (blockName) {
            if (!blockValues) {
                conditionValues.remove(blockName)
            }
            else {
                conditionValues[blockName] = blockValues
            }
        }
        return this
    }

    /**
     * @return
     *      When a conditional block is registered this method returns the collection
     *      of block names visited during the parsing
     */
    Set<String> getConditionalBlockNames() {
        Collections.unmodifiableSet(conditionalNames)
    }

    /**
     * Returns the profile names defined in the config file
     *
     * @return The set of profile names.
     */
    Set<String> getProfileNames() { profileNames }

    private Grengine getGrengine() {
        if( grengine ) {
            return grengine
        }

        // set the required base script
        def config = new CompilerConfiguration()
        config.scriptBaseClass = ConfigBase.class.name
        def params = [:]
        if( renderClosureAsString )
            params.put('renderClosureAsString', true)
        config.addCompilationCustomizers(new ASTTransformationCustomizer(params, ConfigTransform))
        config.addCompilationCustomizers(new ASTTransformationCustomizer(NextflowXform))
        //  add implicit types
        def importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( Duration.name )
        importCustomizer.addImports( MemoryUnit.name )
        config.addCompilationCustomizers(importCustomizer)
        grengine = new Grengine(config)
    }

    ConfigParser setRenderClosureAsString(boolean value) {
        this.renderClosureAsString = value
        return this
    }

    /**
     * Sets any additional variables that should be placed into the binding when evaluating Config scripts
     */
    ConfigParser setBinding(Map vars) {
        this.bindingVars = vars
        return this
    }

    ConfigParser setParams(Map vars) {
        this.paramVars = vars
        return this
    }


    /**
     * Creates a unique name for the config class in order to avoid collision
     * with top level configuration scopes
     *
     * @param text
     * @return
     */
    private String createUniqueName(String text) {
        def hash = Hashing
                .murmur3_32()
                .newHasher()
                .putUnencodedChars(text)
                .hash()
        return "_nf_config_$hash"
    }

    private Script loadScript(String text)  {
        (Script)getGrengine().load(text, createUniqueName(text)).newInstance()
    }

    /**
     * Parses a ConfigObject instances from an instance of java.util.Properties
     * @param The java.util.Properties instance
     */
    ConfigObject parse(Properties properties) {
        ConfigObject config = new ConfigObject()
        for (key in properties.keySet()) {
            def tokens = key.split(/\./)

            def current = config
            def last
            def lastToken
            def foundBase = false
            for (token in tokens) {
                if (foundBase) {
                    // handle not properly nested tokens by ignoring
                    // hierarchy below this point
                    lastToken += "." + token
                    current = last
                } else {
                    last = current
                    lastToken = token
                    current = current."${token}"
                    if (!(current instanceof ConfigObject)) foundBase = true
                }
            }

            if (current instanceof ConfigObject) {
                if (last[lastToken]) {
                    def flattened = last.flatten()
                    last.clear()
                    flattened.each { k2, v2 -> last[k2] = v2 }
                    last[lastToken] = properties.get(key)
                }
                else {
                    last[lastToken] = properties.get(key)
                }
            }
            current = config
        }
        return config
    }
    /**
     * Parse the given script as a string and return the configuration object
     *
     * @see ConfigParser#parse(groovy.lang.Script)
     */
    ConfigObject parse(String text) {
        return parse(loadScript(text))
    }

    /**
     * Parse the given script into a configuration object (a Map)
     * (This method creates a new class to parse the script each time it is called.)
     * @param script The script to parse
     * @return A Map of maps that can be navigating with dot de-referencing syntax to obtain configuration entries
     */
    @Deprecated
    ConfigObject parse(Script script) {
        return parse(script, null)
    }

    /**
     * Parses a Script represented by the given URL into a ConfigObject
     *
     * @param location The location of the script to parse
     * @return The ConfigObject instance
     */
    @Deprecated
    ConfigObject parse(URL location) {
        return parse(loadScript(location.text), FileHelper.asPath(location.toURI()))
    }

    ConfigObject parse(File file) {
        return parse(file.toPath())
    }

    ConfigObject parse(Path path) {
        return parse(loadScript(path.text), path)
    }

    /**
     * Parses the passed groovy.lang.Script instance using the second argument to allow the ConfigObject
     * to retain an reference to the original location other Groovy script
     *
     * @param script The groovy.lang.Script instance
     * @param location The original location of the Script as a URL
     * @return The ConfigObject instance
     */
    ConfigObject parse(Script _script, Path location) {
        final script = (ConfigBase)_script
        Stack<String> currentConditionalBlock = new Stack<String>()
        def config = location ? new ConfigObject(location.toUri().toURL()) : new ConfigObject()
        GroovySystem.metaClassRegistry.removeMetaClass(script.class)
        def mc = script.class.metaClass
        def prefix = ""
        LinkedList stack = new LinkedList()
        LinkedList profileStack = new LinkedList()
        stack << [config: config, scope: [:]]
        def pushStack = { co ->
            stack << [config: co, scope: stack.last.scope.clone()]
        }
        def assignName = { name, co ->
            def current = stack.last
            current.config[name] = co
            current.scope[name] = co
        }
        mc.getProperty = { String name ->
            def current = stack.last
            def result
            if (current.config.get(name)) {
                result = current.config.get(name)
            } else if (current.scope[name]) {
                result = current.scope[name]
            } else {
                try {
                    result = InvokerHelper.getProperty(this, name)
                } catch (GroovyRuntimeException e) {
                    result = new ConfigObject()
                    assignName.call(name, result)
                }
            }
            if( name=='params' && result instanceof Map && paramVars ) {
                result.putAll(paramVars)
            }
            result
        }

        ConfigObject overrides = new ConfigObject()
        mc.invokeMethod = { String name, args ->
            def result
            if (args.length == 1 && args[0] instanceof Closure) {
                if( profileStack && profileStack.last == 'profiles' )
                    profileNames.add(name)

                if (name in conditionValues.keySet()) {
                    try {
                        currentConditionalBlock.push(name)
                        conditionalBlocks.push([:])
                        args[0].call()
                    } finally {
                        currentConditionalBlock.pop()
                        for (entry in conditionalBlocks.pop().entrySet()) {
                            def c = stack.last.config
                            (c != config? c : overrides).merge(entry.value)
                        }
                    }
                } else if (currentConditionalBlock.size() > 0) {
                    String conditionalBlockKey = currentConditionalBlock.peek()
                    conditionalNames.add(name)
                    if (name in conditionValues[conditionalBlockKey]) {
                        def co = conditionalBlocks.peek()[conditionalBlockKey]
                        if( co == null ) {
                            co = new ConfigObject()
                            conditionalBlocks.peek()[conditionalBlockKey] = co
                        }

                        pushStack.call(co)
                        try {
                            currentConditionalBlock.pop()
                            args[0].call()
                        } finally {
                            currentConditionalBlock.push(conditionalBlockKey)
                        }
                        stack.removeLast()
                    }
                } else {
                    def current = stack.last
                    def co
                    if (current.config.get(name) instanceof ConfigObject) {
                        co = current.config.get(name)
                    }
                    else if (current.scope.get(name) instanceof ConfigObject) {
                        co = current.scope.get(name)
                    }
                    else {
                        co = new ConfigObject()
                    }

                    profileStack.add(name)
                    assignName.call(name, co)
                    pushStack.call(co)
                    args[0].call()
                    stack.removeLast()
                    profileStack.removeLast()
                }
            } else if (args.length == 2 && args[1] instanceof Closure) {
                try {
                    prefix = name + '.'
                    assignName.call(name, args[0])
                    args[1].call()
                } finally { prefix = "" }
            } else {
                MetaMethod mm = mc.getMetaMethod(name, args)
                if (mm) {
                    result = mm.invoke(delegate, args)
                } else {
                    throw new MissingMethodException(name, getClass(), args)
                }
            }
            result
        }
        script.metaClass = mc

        def setProperty = { String name, value ->
            assignName.call(prefix + name, value)
        }
        def binding = new ConfigBinding(setProperty)
        if (this.bindingVars) {
            binding.getVariables().putAll(this.bindingVars)
        }

        // add the script file location into the binding
        if( location ) {
            script.setConfigPath(location)
        }

        // disable include parsing when required
        script.setIgnoreIncludes(ignoreIncludes)
        script.setRenderClosureAsString(renderClosureAsString)

        // -- set the binding and run
        script.binding = binding
        script.run()
        config.merge(overrides)

        return config
    }

    /**
     * Disable parsing of {@code includeConfig} directive
     *
     * @param value A boolean value, when {@code true} includes are disabled
     * @return The {@link ConfigParser} object itself
     */
    ConfigParser setIgnoreIncludes(boolean value) {
        this.ignoreIncludes = value
        return this
    }

    /**
     * Since Groovy Script doesn't support overriding setProperty, we have to using a trick with the Binding to provide this
     * functionality
     */
    @PackageScope
    static class ConfigBinding extends Binding {
        Closure callable

        ConfigBinding(Closure c) {
            this.callable = c
        }

        void setVariable(String name, Object value) {
            callable(name, value)
        }
    }
}

