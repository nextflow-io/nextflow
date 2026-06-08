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

package nextflow.cli

import java.lang.reflect.Field

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic

/**
 * Build a machine-readable (JSON) description of the Nextflow command line
 * interface by introspecting the JCommander {@link Parameter} / {@link Parameters}
 * annotations.
 *
 * This powers the global {@code -help-json} flag, which lets tools and LLMs
 * discover the CLI one command at a time without scraping the rendered
 * {@code -help} text:
 *
 * <pre>
 *   nextflow -help-json            # global options + index of all commands
 *   nextflow run -help-json        # full option/argument detail for `run`
 * </pre>
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@CompileStatic
class CliSchema {

    /**
     * Build the schema for the top-level program: the global options plus a
     * name + description index of every available command.
     */
    static String root(CliOptions options, List<CmdBase> commands) {
        final schema = [
            name        : 'nextflow',
            path        : 'nextflow',
            usage       : 'nextflow [options] COMMAND [arg...]',
            params      : paramsOf(CliOptions),
            subcommands : commandIndex(commands),
        ]
        return render(schema)
    }

    /**
     * Build the schema for a single command: its usage, options and positional
     * arguments.
     */
    static String command(CmdBase cmd) {
        final clazz = cmd.getClass()
        final params = paramsOf(clazz)
        final schema = [
            name   : cmd.name,
            path   : "nextflow ${cmd.name}".toString(),
            usage  : usageOf(cmd, params),
            params : params,
        ] as Map<String,Object>

        final description = descriptionOf(clazz)
        if( description )
            schema.help = description
        final aliases = aliasesOf(clazz)
        if( aliases )
            schema.aliases = aliases

        return render(schema)
    }

    /**
     * Map each command name to its description and aliases. Nextflow commands
     * are flat (no nested groups), so a single index lists them all; drill into
     * any one with {@code nextflow <command> -help-json}.
     */
    @CompileDynamic
    private static Map<String,Object> commandIndex(List<CmdBase> commands) {
        final index = new TreeMap<String,Object>()
        for( CmdBase cmd : commands ) {
            final clazz = cmd.getClass()
            final description = descriptionOf(clazz)
            // only advertise commands that opt in with a description, matching `-help`
            if( !description )
                continue
            final entry = [help: description] as Map<String,Object>
            final aliases = aliasesOf(clazz)
            if( aliases )
                entry.aliases = aliases
            index[cmd.name] = entry
        }
        return index
    }

    /**
     * Introspect the {@link Parameter} / {@link DynamicParameter} annotations of
     * a command (or options) class into a list of JSON-friendly maps. Inherited
     * fields are included; the meta options ({@code -h}, {@code -help},
     * {@code -help-json}) are skipped, and hidden options are kept but flagged
     * with {@code hidden: true}.
     */
    @CompileDynamic
    private static List<Map<String,Object>> paramsOf(Class clazz) {
        final instance = newInstance(clazz)
        final result = []
        for( Field field : declaredFields(clazz) ) {
            final p = field.getAnnotation(Parameter)
            final dp = field.getAnnotation(DynamicParameter)
            if( p )
                addParam(result, field, p.names() as List<String>, p.description(), p.hidden(), p.required(), p.arity(), instance)
            else if( dp )
                addParam(result, field, dp.names() as List<String>, dp.description(), dp.hidden(), false, -1, instance)
        }
        // sort for a stable, deterministic ordering: reflection field order is not
        // guaranteed by the JVM, but this is a machine-readable contract
        result.sort { it.opts ? it.opts[0] : it.name }
        return result
    }

    @CompileDynamic
    private static void addParam(List result, Field field, List<String> names, String description, boolean hidden, boolean required, int arity, Object instance) {
        if( names.any { it in META_OPTIONS } )
            return

        final isArgument = !names
        final isFlag = field.getType() == boolean || field.getType() == Boolean || arity == 0
        final entry = [
            name : field.getName(),
            kind : isArgument ? 'argument' : 'option',
            type : typeName(field),
        ] as Map<String,Object>
        if( names )
            entry.opts = names
        if( description )
            entry.help = description
        if( required )
            entry.required = true
        if( isFlag )
            entry.is_flag = true
        // hidden params are kept (parity with rich-click) but flagged so consumers can skip them
        if( hidden )
            entry.hidden = true

        final defValue = defaultOf(field, instance)
        if( defValue != null && !isFlag )
            entry.default = defValue

        result << entry
    }

    /** Build a usage string, e.g. {@code nextflow run [options] <args>}. */
    private static String usageOf(CmdBase cmd, List<Map<String,Object>> params) {
        final pieces = ["nextflow ${cmd.name}".toString()]
        if( params.any { it.kind == 'option' } )
            pieces << '[options]'
        if( params.any { it.kind == 'argument' } )
            pieces << '[args...]'
        return pieces.join(' ')
    }

    @CompileDynamic
    private static String descriptionOf(Class clazz) {
        return clazz.getAnnotation(Parameters)?.commandDescription() ?: null
    }

    @CompileDynamic
    private static List<String> aliasesOf(Class clazz) {
        final names = clazz.getAnnotation(Parameters)?.commandNames()
        return names ? (names as List<String>) : null
    }

    private static String typeName(Field field) {
        return field.getType().getSimpleName()
    }

    /** Read a field's default value from a fresh instance, ignoring empty values. */
    @CompileDynamic
    private static Object defaultOf(Field field, Object instance) {
        if( instance == null )
            return null
        try {
            field.setAccessible(true)
            final value = field.get(instance)
            if( value == null || value == false )
                return null
            if( value instanceof Collection && value.isEmpty() )
                return null
            if( value instanceof Map && value.isEmpty() )
                return null
            return value instanceof CharSequence || value instanceof Number || value instanceof Boolean
                ? value
                : value.toString()
        }
        catch( Exception e ) {
            return null
        }
    }

    private static Object newInstance(Class clazz) {
        try {
            return clazz.getDeclaredConstructor().newInstance()
        }
        catch( Exception e ) {
            return null
        }
    }

    /** Collect declared fields from the class and its superclasses. */
    private static List<Field> declaredFields(Class clazz) {
        final fields = new ArrayList<Field>()
        Class current = clazz
        while( current != null && current != Object ) {
            fields.addAll(current.getDeclaredFields())
            current = current.getSuperclass()
        }
        return fields
    }

    private static String render(Map schema) {
        return JsonOutput.prettyPrint(JsonOutput.toJson(schema))
    }

    private static final List<String> META_OPTIONS = ['-h', '-help', '-help-json']

}
