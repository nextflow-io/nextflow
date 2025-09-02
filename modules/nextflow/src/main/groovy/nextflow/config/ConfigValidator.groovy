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
    private static final List<String> HIDDEN_OPTIONS = List.of(
        'cacheable',
        'dumpChannels',
        'dumpHashes',
        'libDir',
        'poolSize',
        'preview',
        'runName',
        'stubRun',
    );

    /**
     * Core plugin scopes which can only be validated when
     * the plugin is loaded.
     */
    private static final List<String> CORE_PLUGIN_SCOPES = List.of(
        'aws',
        'azure',
        'google',
        'k8s',
        'tower',
        'wave'
    )

    /**
     * Additional config scopes added by third-party plugins
     */
    private SchemaNode.Scope pluginScopes

    ConfigValidator() {
        loadPluginScopes()
    }

    private void loadPluginScopes() {
        final children = new HashMap<String, SchemaNode>()
        for( final scope : Plugins.getExtensions(ConfigScope) ) {
            final clazz = scope.getClass()
            final name = clazz.getAnnotation(ScopeName)?.value()
            final description = clazz.getAnnotation(Description)?.value()
            if( name == '' ) {
                children.putAll(SchemaNode.Scope.of(clazz, '').children())
                continue
            }
            if( !name )
                continue
            if( name in children ) {
                log.warn "Plugin config scope `${clazz.name}` conflicts with existing scope: `${name}`"
                continue
            }
            children.put(name, SchemaNode.Scope.of(clazz, description))
        }
        pluginScopes = new SchemaNode.Scope('', children)
    }

    void validate(ConfigMap config) {
        validate(config.toConfigObject())
    }

    void validate(ConfigObject config) {
        final flatConfig = config.flatten()
        for( String key : flatConfig.keySet() ) {
            final names = key.tokenize('.').findAll { name -> !isSelector(name) }
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
                continue
            if( isValid(names) )
                continue
            if( isMissingCorePluginScope(names.first()) )
                continue
            if( isMapOption(names) )
                continue
            log.warn1 "Unrecognized config option '${fqName}'"
        }
    }

    /**
     * Determine whether a scope name is a process selector.
     *
     * @param name
     */
    private boolean isSelector(String name) {
        return name.startsWith('withLabel:') || name.startsWith('withName:')
    }

    /**
     * Determine whether a config option is defined in the schema.
     *
     * @param names
     */
    boolean isValid(List<String> names) {
        if( names.size() == 1 && names.first() in HIDDEN_OPTIONS )
            return true
        final child = SchemaNode.ROOT.getChild(names)
        if( child instanceof SchemaNode.Option || child instanceof SchemaNode.DslOption )
            return true
        if( pluginScopes.getOption(names) )
            return true
        return false
    }

    /**
     * Determine whether a config scope is from a core plugin
     * which is not currently loaded.
     *
     * @param name
     */
    private boolean isMissingCorePluginScope(String name) {
        return name in CORE_PLUGIN_SCOPES
            && !pluginScopes.children().containsKey(name)
    }

    /**
     * Determine whether a config option is a map option or a
     * property thereof.
     *
     * @param names Config option split into individual names, e.g. 'process.resourceLimits' -> [process, resourceLimits]
     */
    private boolean isMapOption(List<String> names) {
        return isMapOption0(SchemaNode.ROOT, names)
            || isMapOption0(pluginScopes, names)
    }

    private static boolean isMapOption0(SchemaNode.Scope scope, List<String> names) {
        SchemaNode node = scope
        for( final name : names ) {
            if( node instanceof SchemaNode.Scope )
                node = node.children().get(name)
            else if( node instanceof SchemaNode.Placeholder )
                node = node.scope()
            else if( node instanceof SchemaNode.Option )
                return node.type() == Map.class
            else
                return false
        }
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
