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

package nextflow.cli

import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.util.regex.Pattern

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import org.apache.commons.lang3.StringUtils
import org.yaml.snakeyaml.Yaml
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ParamsCollector {

    private Map<String,String> params

    private String paramsFile

    private String stdin

    ParamsCollector(Map<String,String> params, String paramsFile, String stdin) {
        this.params = params
        this.paramsFile = paramsFile ?: SysEnv.get('NXF_PARAMS_FILE')
        this.stdin = stdin
    }

    Map apply(Map configVars) {

        final result = [:]

        // apply params from stdin
        if( stdin ) {
            final text = configVars ? replaceVars0(stdin, configVars) : stdin
            final stdinParams = new JsonSlurper().parseText(text) as Map<String,Object>
            stdinParams.forEach((name, value) -> {
                addParam0(result, name, value)
            })
        }

        // apply params file
        if( paramsFile ) {
            final path = validateParamsFile(paramsFile)
            final type = path.extension.toLowerCase() ?: null
            if( type == 'json' )
                readJsonFile(path, configVars, result)
            else if( type == 'yml' || type == 'yaml' )
                readYamlFile(path, configVars, result)
        }

        // apply CLI params
        for( final entry : params )
            addParam(result, entry.key, entry.value)

        return result
    }

    private static final Pattern DOT_ESCAPED = ~/\\\./
    private static final Pattern DOT_NOT_ESCAPED = ~/(?<!\\)\./

    protected static void addParam(Map params, String key, String value, List path=[], String fullKey=null) {
        if( !fullKey )
            fullKey = key
        final m = DOT_NOT_ESCAPED.matcher(key)
        if( m.find() ) {
            final p = m.start()
            final root = key.substring(0, p)
            if( !root ) throw new AbortOperationException("Invalid parameter name: $fullKey")
            path.add(root)
            def nested = params.get(root)
            if( nested == null ) {
                nested = new LinkedHashMap<>()
                params.put(root, nested)
            }
            else if( nested !instanceof Map ) {
                log.warn "Command line parameter --${path.join('.')} is overwritten by --${fullKey}"
                nested = new LinkedHashMap<>()
                params.put(root, nested)
            }
            addParam((Map)nested, key.substring(p+1), value, path, fullKey)
        }
        else {
            addParam0(params, key.replaceAll(DOT_ESCAPED,'.'), parseParamValue(value))
        }
    }

    protected static void addParam0(Map params, String key, Object value) {
        if( key.contains('-') )
            key = kebabToCamelCase(key)
        params.put(key, value)
    }

    protected static String kebabToCamelCase(String str) {
        final result = new StringBuilder()
        str.split('-').eachWithIndex { String entry, int i ->
            result << (i>0 ? StringUtils.capitalize(entry) : entry )
        }
        return result.toString()
    }

    protected static parseParamValue(String str) {
        if( SysEnv.get('NXF_DISABLE_PARAMS_TYPE_DETECTION') || NF.isSyntaxParserV2() )
            return str

        if( str == null ) return null

        if( str.toLowerCase() == 'true') return Boolean.TRUE
        if( str.toLowerCase() == 'false' ) return Boolean.FALSE

        if( str==~/-?\d+(\.\d+)?/ && str.isInteger() ) return str.toInteger()
        if( str==~/-?\d+(\.\d+)?/ && str.isLong() ) return str.toLong()
        if( str==~/-?\d+(\.\d+)?/ && str.isDouble() ) return str.toDouble()

        return str
    }

    public static final List<String> VALID_PARAMS_FILE = ['json', 'yml', 'yaml']

    private Path validateParamsFile(String file) {
        final result = FileHelper.asPath(file)
        final ext = result.getExtension()
        if( !VALID_PARAMS_FILE.contains(ext) )
            throw new AbortOperationException("Not a valid params file extension: $file -- It must be one of the following: ${VALID_PARAMS_FILE.join(',')}")
        return result
    }

    private static final Pattern PARAMS_VAR = ~/(?m)\$\{(\p{javaJavaIdentifierStart}\p{javaJavaIdentifierPart}*)}/

    protected String replaceVars0(String content, Map binding) {
        content.replaceAll(PARAMS_VAR) { List<String> matcher ->
            // - the regex matcher is represented as list
            // - the first element is the matching string ie. `${something}`
            // - the second element is the group content ie. `something`
            // - make sure the regex contains at least a group otherwise the closure
            // parameter is a string instead of a list of the call fail
            final placeholder = matcher.get(0)
            final key = matcher.get(1)

            if( !binding.containsKey(key) )
                throw new AbortOperationException("Missing params file variable: $placeholder")

            return binding.get(key)
        }
    }

    private void readJsonFile(Path file, Map configVars, Map result) {
        try {
            final text = configVars ? replaceVars0(file.text, configVars) : file.text
            final json = (Map<String,Object>) new JsonSlurper().parseText(text)
            json.forEach((name, value) -> {
                addParam0(result, name, value)
            })
        }
        catch( NoSuchFileException | FileNotFoundException e ) {
            throw new AbortOperationException("Specified params file does not exist: ${file.toUriString()}")
        }
        catch( Exception e ) {
            throw new AbortOperationException("Cannot parse params file: ${file.toUriString()} - Cause: ${e.message}", e)
        }
    }

    private void readYamlFile(Path file, Map configVars, Map result) {
        try {
            final text = configVars ? replaceVars0(file.text, configVars) : file.text
            final yaml = (Map<String,Object>) new Yaml().load(text)
            yaml.forEach((name, value) -> {
                addParam0(result, name, value)
            })
        }
        catch( NoSuchFileException | FileNotFoundException e ) {
            throw new AbortOperationException("Specified params file does not exist: ${file.toUriString()}")
        }
        catch( Exception e ) {
            throw new AbortOperationException("Cannot parse params file: ${file.toUriString()}", e)
        }
    }

}
