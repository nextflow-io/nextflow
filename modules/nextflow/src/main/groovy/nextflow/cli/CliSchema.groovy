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
 *   nextflow -help-json            # global options + the full command tree
 *   nextflow run -help-json        # full option/argument detail for `run`
 * </pre>
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@CompileStatic
class CliSchema {

    /**
     * Build the schema for the top-level program: the global options plus a
     * recursive index of the whole command tree by name + description (e.g.
     * {@code module} → {@code create}, {@code install}, ...). This is the
     * progressive-disclosure entry point — full option/argument detail is
     * revealed by drilling into a command with {@code nextflow <cmd> -help-json}.
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
     * arguments. Commands that dispatch to nested sub-commands (those
     * implementing {@link SubcommandAware}) describe them recursively under a
     * {@code subcommands} key.
     */
    static String command(CmdBase cmd) {
        return render(commandSchema(cmd.name, "nextflow ${cmd.name}".toString(), cmd))
    }

    /**
     * Build the schema map for a command at the given {@code path} (e.g.
     * {@code nextflow module create}). Recurses into nested sub-commands.
     */
    @CompileStatic
    private static Map<String,Object> commandSchema(String name, String path, CmdBase cmd) {
        final clazz = cmd.getClass()
        final params = paramsOf(clazz)
        final schema = [
            name   : name,
            path   : path,
            usage  : usageOf(path, params),
            params : params,
        ] as Map<String,Object>

        final description = descriptionOf(clazz)
        if( description )
            schema.help = description
        final aliases = aliasesOf(clazz)
        if( aliases )
            schema.aliases = aliases
        if( cmd instanceof SubcommandAware ) {
            final subcommands = subcommandIndex(path, (cmd as SubcommandAware).getSubcommands())
            if( subcommands )
                schema.subcommands = subcommands
        }

        return schema
    }

    /**
     * Build the nested {@code subcommands} index for a {@link SubcommandAware}
     * command. Sub-commands backed by a full {@link CmdBase} are introspected
     * recursively for complete option/argument detail; those whose arguments
     * are parsed manually surface a lightweight name + description entry,
     * matching the top-level command index.
     */
    @CompileStatic
    private static Map<String,Object> subcommandIndex(String parentPath, List<SubcommandAware.Subcommand> subcommands) {
        final index = new TreeMap<String,Object>()
        for( SubcommandAware.Subcommand sub : subcommands ) {
            index[sub.name] = sub.command != null
                ? commandSchema(sub.name, "${parentPath} ${sub.name}".toString(), sub.command)
                : leafEntry(sub.help, sub.aliases)
        }
        return index
    }

    /** A lightweight index entry carrying only help text and aliases, when present. */
    @CompileDynamic
    private static Map<String,Object> leafEntry(String help, List<String> aliases) {
        final entry = [:] as Map<String,Object>
        if( help )
            entry.help = help
        if( aliases )
            entry.aliases = aliases
        return entry
    }

    /**
     * Lightweight, recursive index of the command tree: each command mapped to
     * its description plus — for commands that dispatch sub-commands — a nested
     * index of the same shape, all the way down. Deliberately carries only names
     * and help text (no params, aliases or paths): it's the map of the CLI, and
     * full detail is disclosed by drilling into a command.
     */
    @CompileDynamic
    private static Map<String,Object> commandIndex(List<CmdBase> commands) {
        final index = new TreeMap<String,Object>()
        for( CmdBase cmd : commands ) {
            final description = descriptionOf(cmd.getClass())
            // only advertise commands that opt in with a description, matching `-help`
            if( !description )
                continue
            index[cmd.name] = indexEntry(description, cmd)
        }
        return index
    }

    /** One node of the lightweight index: help text, plus a nested index when the command dispatches sub-commands. */
    @CompileDynamic
    private static Map<String,Object> indexEntry(String help, CmdBase cmd) {
        final entry = [:] as Map<String,Object>
        if( help )
            entry.help = help
        if( cmd instanceof SubcommandAware ) {
            final nested = new TreeMap<String,Object>()
            for( SubcommandAware.Subcommand sub : (cmd as SubcommandAware).getSubcommands() ) {
                final subHelp = sub.command != null ? descriptionOf(sub.command.getClass()) : sub.help
                nested[sub.name] = sub.command != null ? indexEntry(subHelp, sub.command) : leafHelp(subHelp)
            }
            if( nested )
                entry.subcommands = nested
        }
        return entry
    }

    /** A leaf index node carrying only help text (for manually-parsed sub-commands). */
    @CompileDynamic
    private static Map<String,Object> leafHelp(String help) {
        final entry = [:] as Map<String,Object>
        if( help )
            entry.help = help
        return entry
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

    /** Build a usage string, e.g. {@code nextflow run [options] [args...]}. */
    private static String usageOf(String path, List<Map<String,Object>> params) {
        final pieces = [path]
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

    /** Soft line-length budget for the compact renderer. */
    private static final int MAX_WIDTH = 100
    private static final int INDENT = 2

    private static String render(Map schema) {
        return format(schema, 0, 0)
    }

    /**
     * Render JSON that stays human-readable but is condensed: any object or
     * array that fits within {@link #MAX_WIDTH} columns is kept on a single
     * line, and only the containers that overflow are expanded — so leaf maps
     * like a single parameter print as one tidy line instead of a dozen.
     *
     * @param value the value to render
     * @param indent leading-whitespace width of the line this value starts on
     * @param col    column at which this value's first character is placed
     *               (accounts for any {@code "key": } prefix already emitted)
     */
    @CompileDynamic
    private static String format(Object value, int indent, int col) {
        final oneLine = inline(value)
        if( col + oneLine.length() <= MAX_WIDTH )
            return oneLine

        final childIndent = indent + INDENT
        final pad = ' ' * childIndent
        if( value instanceof Map && !value.isEmpty() ) {
            final parts = value.collect { k, v ->
                final key = JsonOutput.toJson(k.toString()) + ': '
                pad + key + format(v, childIndent, childIndent + key.length())
            }
            return '{\n' + parts.join(',\n') + '\n' + (' ' * indent) + '}'
        }
        if( value instanceof List && !value.isEmpty() ) {
            final parts = value.collect { v -> pad + format(v, childIndent, childIndent) }
            return '[\n' + parts.join(',\n') + '\n' + (' ' * indent) + ']'
        }
        // scalar (or empty container) that can't be broken — emit as-is
        return oneLine
    }

    /** Single-line JSON for a value, with {@code ": "} / {@code ", "} spacing for readability (unlike minified output). */
    @CompileDynamic
    private static String inline(Object value) {
        if( value instanceof Map )
            return value.isEmpty() ? '{}' : '{' + value.collect { k, v -> JsonOutput.toJson(k.toString()) + ': ' + inline(v) }.join(', ') + '}'
        if( value instanceof List )
            return value.isEmpty() ? '[]' : '[' + value.collect { inline(it) }.join(', ') + ']'
        return JsonOutput.toJson(value)
    }

    private static final List<String> META_OPTIONS = ['-h', '-help', '-help-json']

}
