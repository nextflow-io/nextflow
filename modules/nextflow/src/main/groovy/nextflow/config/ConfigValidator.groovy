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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.config.schema.ConfigScope
import nextflow.config.schema.SchemaNode
import nextflow.config.schema.ScopeName
import nextflow.plugin.Plugins
import nextflow.script.dsl.Description
/**
 * Validate the Nextflow configuration
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ConfigValidator {

    /**
     * Hidden options added by ConfigBuilder
     */
    private static final List<String> hiddenOptions = List.of(
        'bucketDir',
        'cacheable',
        'dumpChannels',
        'libDir',
        'poolSize',
        'preview',
        'runName',
        'stubRun',
    );

    /**
     * Additional config scopes added by third-party plugins
     */
    private SchemaNode.Scope pluginScopes

    ConfigValidator() {
        loadPluginScopes()
    }

    private void loadPluginScopes() {
        final scopes = new HashMap<String, SchemaNode>()
        for( final scope : Plugins.getExtensions(ConfigScope) ) {
            final clazz = scope.getClass()
            final name = clazz.getAnnotation(ScopeName)?.value()
            final description = clazz.getAnnotation(Description)?.value()
            if( !name )
                continue
            if( name in scopes ) {
                log.warn "Plugin config scope `${clazz.name}` conflicts with existing scope: `${name}`"
                continue
            }
            scopes.put(name, SchemaNode.Scope.of(clazz, description))
        }
        pluginScopes = new SchemaNode.Scope('', scopes)
    }

    void validate(ConfigMap config) {
        if( NF.getSyntaxParserVersion() != 'v2' )
            return
        validate(config.toConfigObject())
    }

    void validate(ConfigObject config) {
        if( NF.getSyntaxParserVersion() != 'v2' )
            return
        final flatConfig = config.flatten()
        for( String key : flatConfig.keySet() ) {
            final names = key.tokenize('.')
            if( names.first() == 'profiles' ) {
                if( !names.isEmpty() ) names.remove(0)
                if( !names.isEmpty() ) names.remove(0)
            }
            final scope = names.first()
            if( scope == 'env' ) {
                checkEnv(names.last())
                continue
            }
            if( scope == 'params' )
                continue
            final fqName = names.join('.')
            if( fqName.startsWith('process.ext.') )
                return
            if( !isValid(names) ) {
                log.warn "Unrecognized config option '${fqName}'"
                continue
            }
        }
    }

    /**
     * Determine whether a config option is defined in the schema.
     *
     * @param names
     */
    boolean isValid(List<String> names) {
        if( names.size() == 1 && names.first() in hiddenOptions )
            return true
        final child = SchemaNode.ROOT.getChild(names)
        if( child instanceof SchemaNode.Option || child instanceof SchemaNode.DslOption )
            return true
        if( pluginScopes.getOption(names) )
            return true
        return false
    }

    /**
     * Warn about setting `NXF_*` environment variables in the config.
     *
     * @param name
     */
    private void checkEnv(String name) {
        if( name.startsWith('NXF_') && name!='NXF_DEBUG' )
            log.warn "Nextflow environment variables must be defined in the launch environment -- the following environment variable in the config will be ignored: '$name'"
    }
}
