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

    private List<String> profiles

    private Map target = [:]

    private Set<String> parsedProfiles = []

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

    void setProfiles(List<String> profiles) {
        this.profiles = profiles
    }

    void addParsedProfile(String profile) {
        parsedProfiles.add(profile)
    }

    Set<String> getParsedProfiles() {
        return parsedProfiles
    }

    Map getTarget() {
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
        final dsl = blockDsl(names)
        final cl = (Closure)closure.clone()
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.setDelegate(dsl)
        cl.call()
        dsl.apply()
    }

    private ConfigBlockDsl blockDsl(List<String> names) {
        if( names.size() == 1 && names.first() == 'plugins' )
            return new PluginsDsl(this)

        if( names.size() == 1 && names.first() == 'process' )
            return new ProcessDsl(this, names)

        if( names.size() == 1 && names.first() == 'profiles' )
            return new ProfilesDsl(this, profiles)

        return new ConfigBlockDsl(this, names)
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
        final parser = new ConfigParserV2()
                .setIgnoreIncludes(ignoreIncludes)
                .setRenderClosureAsString(renderClosureAsString)
                .setStrict(strict)
                .setBinding(binding.getVariables())
                .setParams(target.params as Map)
                .setProfiles(profiles)
        final config = parser.parse(configText, includePath)
        parsedProfiles.addAll(parser.getProfiles())

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

    private static class ConfigBlockDsl {
        protected ConfigDsl dsl
        protected List<String> scope

        ConfigBlockDsl(ConfigDsl dsl, List<String> scope) {
            this.dsl = dsl
            this.scope = scope
        }

        void assign(List<String> names, Object right) {
            dsl.assign(scope + names, right)
        }

        void block(String name, Closure closure) {
            dsl.block(scope + [name], closure)
        }

        void withLabel(String label, Closure closure) {
            throw new ConfigParseException("Process selectors are only allowed in the `process` scope (offending scope: `${scope.join('.')}`)")
        }

        void withName(String selector, Closure closure) {
            throw new ConfigParseException("Process selectors are only allowed in the `process` scope (offending scope: `${scope.join('.')}`)")
        }

        void includeConfig(String includeFile) {
            dsl.includeConfig(scope, includeFile)
        }

        void apply() {
        }
    }

    private static class PluginsDsl extends ConfigBlockDsl {
        PluginsDsl(ConfigDsl dsl) {
            super(dsl, Collections.emptyList())
        }

        void id(String value) {
            final target = dsl.getTarget()
            final plugins = (Set) target.computeIfAbsent('plugins', (k) -> new HashSet<>())
            plugins.add(value)
        }
    }

    private static class ProcessDsl extends ConfigBlockDsl {
        ProcessDsl(ConfigDsl dsl, List<String> scope) {
            super(dsl, scope)
        }

        @Override
        void withLabel(String label, Closure closure) {
            dsl.block(scope + ["withLabel:${label}".toString()], closure)
        }

        @Override
        void withName(String selector, Closure closure) {
            dsl.block(scope + ["withName:${selector}".toString()], closure)
        }
    }

    private static class ProfilesDsl extends ConfigBlockDsl {
        private List<String> profiles
        private Map<String,Closure> blocks = [:]

        ProfilesDsl(ConfigDsl dsl, List<String> profiles) {
            super(dsl, Collections.emptyList())
            this.profiles = profiles
        }

        @Override
        void assign(List<String> names, Object right) {
            throw new ConfigParseException("Only profile blocks are allowed in the `profiles` scope")
        }

        @Override
        void block(String name, Closure closure) {
            blocks[name] = closure
            dsl.addParsedProfile(name)
        }

        @Override
        void includeConfig(String includeFile) {
            throw new ConfigParseException("Only profile blocks are allowed in the `profiles` scope")
        }

        @Override
        void apply() {
            if( profiles != null ) {
                // apply profiles in the order they were specified
                for( final name : profiles ) {
                    final closure = blocks[name]
                    if( closure )
                        dsl.block(scope, closure)
                }
            }
            else {
                // append all profiles to the config map
                for( final name : blocks.keySet() ) {
                    final closure = blocks[name]
                    dsl.block(['profiles', name], closure)
                }
            }
        }
    }

}
