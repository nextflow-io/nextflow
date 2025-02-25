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
import nextflow.config.dsl.ConfigSchema
import nextflow.config.dsl.ConfigScope
import nextflow.plugin.Plugins
/**
 * Validate the Nextflow configuration
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ConfigValidator {

    void validate(ConfigMap config) {
        validate(config.toConfigObject())
    }

    void validate(ConfigObject config) {
        final schema = getSchema()
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
            if( fqName !in schema ) {
                log.warn "Unrecognized config option '${fqName}'"
                continue
            }
        }
    }

    protected Set<String> getSchema() {
        final schema = new HashSet<String>()
        schema.addAll(ConfigSchema.OPTIONS.keySet())
        for( final scope : Plugins.getExtensions(ConfigScope) )
            schema.addAll(ConfigSchema.getConfigOptions(scope).keySet())
        // hidden options added by ConfigBuilder
        schema.addAll(List.of(
            'bucketDir',
            'cacheable',
            'dumpChannels',
            'libDir',
            'poolSize',
            'plugins',
            'preview',
            'runName',
            'stubRun',
        ))
        return schema
    }

    /**
     * Warn about setting `NXF_*` environment variables in the config.
     *
     * @param name
     */
    protected void checkEnv(String name) {
        if( name.startsWith('NXF_') && name!='NXF_DEBUG' )
            log.warn "Nextflow environment variables must be defined in the launch environment -- the following environment variable in the config will be ignored: `$name`"
    }
}