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

package nextflow.config.parser.v2

import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.exception.ConfigParseException
import nextflow.extension.Bolts
import nextflow.file.FileHelper
/**
 * Builder DSL for Nextflow config files.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ConfigDsl extends Script {

    private boolean ignoreIncludes

    private boolean renderClosureAsString

    private boolean strict

    private Path configPath

    private Map target = [:]

    void setIgnoreIncludes(boolean value) {
        this.ignoreIncludes = value
    }

    void setRenderClosureAsString(boolean value) {
        this.renderClosureAsString = value
    }

    void setStrict(boolean value) {
        this.strict = value
    }

    void setConfigPath(Path path) {
        this.configPath = path
    }

    void setParams(Map params) {
        target.params = params
    }

    Map getTarget() {
        if( !target.params )
            target.remove('params')
        return target
    }

    Object run() {}

    @Override
    def getProperty(String name) {
        if( name == 'params' )
            return target.params

        try {
            return super.getProperty(name)
        }
        catch( MissingPropertyException e ) {
            if( strict )
                throw e
            else
                return null
        }
    }

    void append(List<String> names, Object right) {
        final values = (Set) navigate(names.init()).computeIfAbsent(names.last(), (k) -> new HashSet<>())
        values.add(right)
    }

    void assign(List<String> names, Object right) {
        navigate(names.init()).put(names.last(), right)
    }

    private Map navigate(List<String> names) {
        Map ctx = target
        for( final name : names ) {
            if( name !in ctx ) ctx[name] = [:]
            ctx = ctx[name] as Map
        }
        return ctx
    }

    void block(String name, Closure closure) {
        block([name], closure)
    }

    void block(List<String> names, Closure closure) {
        final delegate = new ConfigBlockDsl(this, names)
        final cl = (Closure)closure.clone()
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.setDelegate(delegate)
        cl.call()
    }

    /**
     * Get the value of an environment variable from the launch environment.
     *
     * @param name
     */
    String env(String name) {
        return SysEnv.get(name)
    }

    void includeConfig(String includeFile) {
        includeConfig([], includeFile)
    }

    void includeConfig(List<String> names, String includeFile) {
        assert includeFile

        if( ignoreIncludes )
            return

        Path includePath = FileHelper.asPath(includeFile)
        log.trace "Include config file: $includeFile [parent: $configPath]"

        if( !includePath.isAbsolute() && configPath )
            includePath = configPath.resolveSibling(includeFile)

        final configText = readConfigFile(includePath)
        final config = new ConfigParserV2()
                .setIgnoreIncludes(ignoreIncludes)
                .setRenderClosureAsString(renderClosureAsString)
                .setStrict(strict)
                .setBinding(binding.getVariables())
                .parse(configText, includePath)

        final ctx = navigate(names)
        ctx.putAll(Bolts.deepMerge(ctx, config))
    }

    /**
     * Read the content of a config file. The result is cached to
     * avoid multiple reads.
     *
     * @param includePath
     */
    @Memoized
    protected static String readConfigFile(Path includePath) {
        try {
            return includePath.getText()
        }
        catch (NoSuchFileException | FileNotFoundException ignored) {
            throw new NoSuchFileException("Config file does not exist: ${includePath.toUriString()}")
        }
        catch (IOException e) {
            throw new IOException("Cannot read config file include: ${includePath.toUriString()}", e)
        }
    }

    static class ConfigBlockDsl {
        private ConfigDsl dsl
        private List<String> scope

        ConfigBlockDsl(ConfigDsl dsl, List<String> scope) {
            this.dsl = dsl
            this.scope = scope
        }

        void append(String name, Object right) {
            dsl.append(scope, right)
        }

        void assign(List<String> names, Object right) {
            dsl.assign(scope + names, right)
        }

        void block(String name, Closure closure) {
            dsl.block(scope + [name], closure)
        }

        void withLabel(String label, Closure closure) {
            if( !isWithinProcessScope() )
                throw new ConfigParseException("Process selectors are only allowed in the `process` scope (offending scope: `${scope.join('.')}`)")
            dsl.block(scope + ["withLabel:${label}".toString()], closure)
        }

        void withName(String selector, Closure closure) {
            if( !isWithinProcessScope() )
                throw new ConfigParseException("Process selectors are only allowed in the `process` scope (offending scope: `${scope.join('.')}`)")
            dsl.block(scope + ["withName:${selector}".toString()], closure)
        }

        private boolean isWithinProcessScope() {
            if( scope.size() == 1 )
                return scope.first() == 'process'
            if( scope.size() == 3 )
                return scope.first() == 'profiles' && scope.last() == 'process'
            return false
        }

        void includeConfig(String includeFile) {
            dsl.includeConfig(scope, includeFile)
        }
    }

}
