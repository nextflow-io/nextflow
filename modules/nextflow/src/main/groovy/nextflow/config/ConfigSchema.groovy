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
import nextflow.plugin.Plugins
import org.pf4j.ExtensionPoint

/**
 * Implements the validation of Nextflow configuration against
 * typed schemas.
 */
@Slf4j
@CompileStatic
class ConfigSchemaValidator {

    private static Map<String,Class> schema = null

    static void validate(ConfigObject config) {
        if( !schema ) {
            schema = [:]
            for( def clazz : Plugins.getExtensions(ConfigSchema) ) {
                for( def field : clazz.getClass().getDeclaredFields() ) {
                    final annot = field.getAnnotation(ConfigOption)
                    if( !annot )
                        continue
                    if( annot.value() in schema )
                        log.warn "Config option '${annot.value()}' is defined multiple times"
                    schema.put(annot.value(), field.getType())

                    log.info "config option '${annot.value()}' -> ${field.getType()}"
                }
                for( def method : clazz.getClass().getDeclaredMethods() ) {
                    final annot = method.getAnnotation(ConfigOption)
                    if( !annot )
                        continue
                    if( annot.value() in schema )
                        log.warn "Config option '${annot.value()}' is defined multiple times"
                    schema.put(annot.value(), method.getReturnType())

                    log.info "config option '${annot.value()}' -> ${method.getReturnType().getSimpleName()}"
                }
            }
        }

        final flatConfig = config.flatten()
        for( String key : flatConfig.keySet() ) {
            final scope = key.tokenize('.').first()
            if( scope == 'params' )
                continue
            if( scope == 'profiles' )
                key = key.tokenize('.')[2..-1].join('.')
            if( key !in schema ) {
                log.warn "Unrecognized config option '${key}'"
                continue
            }
            final value = flatConfig[key]
            if( value == null )
                continue
            final expectedType = schema[key]
            final actualType = value.getClass()
            if( !expectedType.isAssignableFrom(actualType) ) {
                log.warn "Invalid config option '${key}' -- expected a ${expectedType} but got a ${actualType}"
                continue
            }
        }
    }

}

/**
 * Interface for classes that define config options
 */
interface ConfigSchema extends ExtensionPoint {}
