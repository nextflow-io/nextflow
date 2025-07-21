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

package nextflow.config

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.exception.AbortOperationException
import nextflow.exception.ConfigParseException
import nextflow.secret.SecretsLoader
/**
 * Builds up the Nextflow configuration object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ConfigBuilder {

    public static final String DEFAULT_PROFILE = 'standard'

    protected Path baseDir

    protected Path currentDir

    protected Map params

    protected String profile = DEFAULT_PROFILE

    protected boolean validateProfile

    protected boolean ignoreIncludes

    protected boolean showAllProfiles

    protected boolean showClosures

    protected boolean showMissingVariables

    protected boolean stripSecrets

    protected List<Path> parsedConfigFiles = []

    protected Map<ConfigObject, String> emptyVariables = [:]

    protected List<String> warnings = []

    ConfigBuilder setBaseDir(Path path) {
        this.baseDir = path.complete()
        return this
    }

    ConfigBuilder setCurrentDir(Path path) {
        this.currentDir = path.complete()
        return this
    }

    ConfigBuilder setIgnoreIncludes(boolean value) {
        this.ignoreIncludes = value
        return this
    }

    ConfigBuilder setParams(Map params) {
        this.params = params
        return this
    }

    ConfigBuilder setProfile(String value) {
        this.profile = value ?: DEFAULT_PROFILE
        this.validateProfile = value as boolean
        return this
    }

    ConfigBuilder setShowAllProfiles(boolean value) {
        this.showAllProfiles = value
        return this
    }

    ConfigBuilder setShowClosures(boolean value) {
        this.showClosures = value
        return this
    }

    ConfigBuilder setShowMissingVariables(boolean value) {
        this.showMissingVariables = value
        return this
    }

    ConfigBuilder setStripSecrets(boolean value) {
        this.stripSecrets = value
        return this
    }

    List<Path> getParsedConfigFiles() {
        return parsedConfigFiles
    }

    List<String> getWarnings() {
        return warnings
    }

    Map configVars() {
        // this is needed to make sure to reuse the same
        // instance of the config vars across different instances of the ConfigBuilder
        // and prevent multiple parsing of the same params file (which can even be remote resource)
        return cacheableConfigVars(baseDir)
    }

    @Memoized
    private static Map cacheableConfigVars(Path base) {
        final binding = new HashMap(10)
        binding.put('baseDir', base)
        binding.put('projectDir', base)
        binding.put('launchDir', Path.of('.').toRealPath())
        binding.put('outputDir', Path.of('results').complete())
        binding.put('secrets', SecretsLoader.secretContext())
        return binding
    }

    /**
     * Build a config object from the given config entries and
     * environment variables. Each entry can be either a file (Path)
     * or a snippet (String).
     *
     * @param env
     * @param configEntries
     */
    ConfigObject build(Map<String,String> env=[:], List configEntries)  {
        assert env != null

        final parser = ConfigParserFactory.create()
                .setRenderClosureAsString(showClosures)
                .setStripSecrets(stripSecrets)
                .setIgnoreIncludes(ignoreIncludes)

        if( params )
            parser.setParams(params)

        final result = new ConfigObject()
        result.put('env', new ConfigObject())

        // add the user specified environment to the session env
        final resultEnv = (Map) result.env
        env.sort().forEach((name, value) -> {
            resultEnv.put(name, value)
        })

        if( configEntries ) {
            final binding = [:]
            // the configuration object binds always the current environment
            // so that in the configuration file may be referenced any variable
            // in the current environment
            if( NF.getSyntaxParserVersion() == 'v1' ) {
                binding.putAll(System.getenv())
                binding.putAll(env)
            }
            binding.putAll(configVars())

            parser.setBinding(binding)

            // merge of the provided configuration files
            for( final entry : configEntries ) {
                try {
                    merge0(result, parser, entry)
                }
                catch( ConfigParseException e ) {
                    throw e
                }
                catch( Exception e ) {
                    final message = entry instanceof Path
                        ? "Unable to parse config file: '$entry'".toString()
                        : "Unable to parse configuration "
                    throw new ConfigParseException(message, e)
                }
            }

            // validate profiles if specified
            if( validateProfile ) {
                checkValidProfile(parser.getProfiles())
            }
        }

        // guarantee top scopes
        for( final name : List.of('env','executor','params','process') ) {
            if( !result.isSet(name) )
                result.put(name, new ConfigObject())
        }

        return result
    }

    /**
     * Merge a config entry into an accumulated config object.
     *
     * @param result The main {@link ConfigObject}
     * @param parser The {@ConfigParser} instance
     * @param entry The next config file or snippet to parse
     */
    protected void merge0(ConfigObject result, ConfigParser parser, entry) {
        if( !entry )
            return

        // select the profile
        if( !showAllProfiles ) {
            log.debug "Applying config profile: `${profile}`"
            parser.setProfiles(profile.tokenize(','))
        }

        final config = parse0(parser,entry)
        if( NF.getSyntaxParserVersion() == 'v1' )
            checkUnresolvedConfig(config,entry)
        result.merge(config)
    }

    /**
     * Parse a config file or snippet.
     *
     * @param parser
     * @param entry
     */
    protected ConfigObject parse0(ConfigParser parser, entry) {
        if( entry instanceof File ) {
            return parse0(parser, entry.toPath())
        }

        if( entry instanceof Path ) {
            parsedConfigFiles << entry
            return parser.parse(entry)
        }

        if( entry instanceof CharSequence ) {
            return parser.parse(entry.toString())
        }

        throw new IllegalStateException("Unexpected config entry: ${entry}")
    }

    /**
     * Verify that a config object does not contain any unresolved attributes.
     *
     * @param config The {@link ConfigObject} to verify
     * @param file The source config file/snippet
     */
    protected void checkUnresolvedConfig(ConfigObject config, file, String parent=null, List<ConfigObject> stack = []) {
        for( String key : new ArrayList<>(config.keySet()) ) {
            final value = config.get(key)
            if( value instanceof ConfigObject ) {
                final fqKey = parent
                    ? "${parent}.${key}".toString()
                    : key as String
                if( value.isEmpty() ) {
                    final msg = "Unknown config attribute `$fqKey` -- check config file: $file".toString()
                    if( showMissingVariables ) {
                        emptyVariables.put(value, key)
                        warnings.add(msg)
                    }
                    else {
                        log.debug("In the following config snippet the attribute `$fqKey` is empty:\n${->config.prettyPrint().indent('  ')}")
                        throw new ConfigParseException(msg)
                    }
                }
                else {
                    stack.push(config)
                    try {
                        if( !stack.contains(value)) {
                            checkUnresolvedConfig(value, file, fqKey, stack)
                        }
                        else {
                            log.debug("Found a recursive config property: `$fqKey`")
                        }
                    }
                    finally {
                        stack.pop()
                    }
                }
            }
            else if( value instanceof GString && showMissingVariables ) {
                final str = (GString) value
                for( int i=0; i<str.values.length; i++ ) {
                    // try replace empty interpolated strings with variable handle
                    final arg = str.values[i]
                    final name = emptyVariables.get(arg)
                    if( name )
                        str.values[i] = '$' + name
                }
            }
        }
    }

    /**
     * Validate that each specified profile exists in the merged config.
     *
     * @param declaredProfiles
     */
    protected void checkValidProfile(Collection<String> declaredProfiles) {
        if( !profile || profile == DEFAULT_PROFILE ) {
            return
        }

        log.debug "Available config profiles: $declaredProfiles"
        for( String name : profile.tokenize(',') ) {
            if( name in declaredProfiles )
                continue

            def message = "Unknown configuration profile: '${name}'"
            def choices = declaredProfiles.closest(name)
            if( choices ) {
                message += "\n\nDid you mean one of these?\n"
                choices.each { message += "    ${it}\n" }
                message += '\n'
            }

            throw new AbortOperationException(message)
        }
    }
}
