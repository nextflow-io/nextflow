/*
 * Copyright 2013-2026, Seqera Labs
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

import java.nio.file.Files
import java.nio.file.Path

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.SysEnv

/**
 * Coerces CLI param string values to typed values using a
 * {@code nextflow_schema.json} file (JSON Schema, as used by nf-core
 * pipelines) as a fallback type source.
 *
 * Under syntax parser v2, CLI params arrive as strings. When the pipeline
 * has not declared typed params in main.nf or non-null defaults in
 * nextflow.config, those strings stay strings -- and may break param logic
 * that expects e.g. numeric comparison. If a {@code nextflow_schema.json}
 * lives next to main.nf, this helper reads property types from it and
 * coerces CLI values accordingly, giving the pipeline typed params for
 * free without requiring main.nf changes.
 *
 * Coercion is best-effort and non-destructive: values that don't match a
 * declared type are left as strings, and a missing or malformed schema
 * is silently ignored.
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@Slf4j
@CompileStatic
class SchemaParamsHelper {

    static final String SCHEMA_FILENAME = 'nextflow_schema.json'

    /**
     * Apply schema-based type coercion in place on a CLI params map.
     *
     * @param baseDir   pipeline project base directory (where main.nf lives)
     * @param cliParams CLI params map; mutated in place
     */
    static void applySchemaTypes(Path baseDir, Map<String,?> cliParams) {
        if( !cliParams || baseDir == null )
            return
        if( !NF.isSyntaxParserV2() )
            return
        if( SysEnv.get('NXF_DISABLE_PARAMS_TYPE_DETECTION') )
            return

        final schemaFile = baseDir.resolve(SCHEMA_FILENAME)
        if( !Files.exists(schemaFile) )
            return

        final types = readSchemaTypes(schemaFile)
        if( !types )
            return

        log.debug "Applying types from ${schemaFile} to ${cliParams.size()} CLI param(s) -- ${types.size()} param type(s) declared in schema"
        coerceInPlace(cliParams, types)
    }

    /**
     * Parse a JSON schema file and return a map of {@code paramName -> jsonType}
     * (e.g. {@code "integer"}, {@code "number"}, {@code "boolean"}).
     */
    static Map<String,String> readSchemaTypes(Path schemaFile) {
        try {
            final root = new JsonSlurper().parse(schemaFile)
            final types = new LinkedHashMap<String,String>()
            collectProperties(root, types)
            return types
        }
        catch( Exception e ) {
            log.warn "Unable to parse ${schemaFile} for fallback param typing -- ${e.message}"
            return Collections.<String,String>emptyMap()
        }
    }

    /**
     * Recursively walk a JSON Schema fragment, collecting top-level property
     * names and their declared {@code type}. Handles nf-core-style schemas
     * that nest properties under {@code definitions} or {@code $defs}, plus
     * any {@code allOf}/{@code oneOf}/{@code anyOf} compositions.
     */
    private static void collectProperties(Object node, Map<String,String> types) {
        if( node !instanceof Map )
            return
        final map = (Map) node

        final props = map.get('properties')
        if( props instanceof Map ) {
            for( final entry in (Map<String,Object>) props ) {
                final schemaMap = entry.value instanceof Map ? (Map) entry.value : null
                if( schemaMap == null )
                    continue
                final type = schemaMap.get('type')
                if( type instanceof String && !types.containsKey(entry.key) )
                    types.put(entry.key, (String) type)
            }
        }

        for( final key in ['definitions', '$defs', 'allOf', 'oneOf', 'anyOf'] ) {
            final sub = map.get(key)
            final children = sub instanceof Map ? ((Map) sub).values()
                : sub instanceof List ? (List) sub
                : null
            if( children == null )
                continue
            for( final child in children )
                collectProperties(child, types)
        }
    }

    private static void coerceInPlace(Map<String,?> params, Map<String,String> types) {
        for( final name : new ArrayList<String>(params.keySet()) ) {
            final value = params.get(name)
            final coerced = coerceValue(value, types.get(name))
            // coerceValue returns the same reference when no coercion applied
            if( coerced !== value )
                ((Map) params).put(name, coerced)
        }
    }

    private static Object coerceValue(Object value, String type) {
        if( value !instanceof CharSequence )
            return value
        if( !type )
            return value
        final str = value.toString()
        switch( type ) {
            case 'boolean':
                if( str.equalsIgnoreCase('true') ) return Boolean.TRUE
                if( str.equalsIgnoreCase('false') ) return Boolean.FALSE
                break
            case 'integer':
                if( str.isInteger() ) return str.toInteger()
                if( str.isLong() ) return str.toLong()
                if( str.isBigInteger() ) return str.toBigInteger()
                break
            case 'number':
                if( str.isInteger() ) return str.toInteger()
                if( str.isLong() ) return str.toLong()
                if( str.isFloat() ) return str.toFloat()
                if( str.isDouble() ) return str.toDouble()
                if( str.isBigDecimal() ) return str.toBigDecimal()
                break
        }
        return value
    }
}
