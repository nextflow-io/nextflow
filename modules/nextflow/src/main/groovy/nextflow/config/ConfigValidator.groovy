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
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.SpecNode
import nextflow.config.spec.ScopeName
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
    private SpecNode.Scope pluginScopes

    ConfigValidator() {
        loadPluginScopes()
    }

    private void loadPluginScopes() {
        final children = new HashMap<String, SpecNode>()
        for( final scope : Plugins.getExtensions(ConfigScope) ) {
            final clazz = scope.getClass()
            final name = clazz.getAnnotation(ScopeName)?.value()
            final description = clazz.getAnnotation(Description)?.value()
            if( name == '' ) {
                children.putAll(SpecNode.Scope.of(clazz, '').children())
                continue
            }
            if( !name )
                continue
            if( name in children ) {
                log.warn "Plugin config scope `${clazz.name}` conflicts with existing scope: `${name}`"
                continue
            }
            children.put(name, SpecNode.Scope.of(clazz, description))
        }
        pluginScopes = new SpecNode.Scope('', children)
    }

    /**
     * Validate a config block within the given scope.
     *
     * @param config
     * @param scopes
     */
    void validate(Map<String,?> config, List<String> scopes=[]) {
        for( final entry : config.entrySet() ) {
            final key = entry.key
            final value = entry.value

            final names = scopes + [key]

            if( names.size() == 2 && names.first() == 'profiles' )
                names.clear()

            if( value instanceof Map ) {
                if( isSelector(key) )
                    names.removeLast()
                if( isMapOption(names) )
                    continue
                validate(value, names)
            }
            else {
                validateOption(names)
            }
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
     * Validate a config option given by the list of names.
     *
     * For example, the option 'process.resourceLimits' is represented
     * as ['process', 'resourceLimits'].
     *
     * @param names
     */
    void validateOption(List<String> names) {
        final scope = names.first()
        if( scope == 'env' ) {
            checkEnv(names.last())
            return
        }
        if( scope == 'params' )
            return
        if( isMissingCorePluginScope(scope) )
            return
        if( isValid(names) )
            return
        log.warn1 "Unrecognized config option '${names.join('.')}'"
    }

    /**
     * Determine whether a config option is defined in the spec.
     *
     * @param names
     */
    boolean isValid(List<String> names) {
        if( names.size() == 1 && names.first() in HIDDEN_OPTIONS )
            return true
        final child = SpecNode.ROOT.getChild(names)
        if( child instanceof SpecNode.Option || child instanceof SpecNode.DslOption )
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
     * @param names
     */
    private boolean isMapOption(List<String> names) {
        return isMapOption0(SpecNode.ROOT, names)
            || isMapOption0(pluginScopes, names)
    }

    private static boolean isMapOption0(SpecNode.Scope scope, List<String> names) {
        final node = scope.getOption(names)
        return node != null && node.types().contains(Map.class)
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
