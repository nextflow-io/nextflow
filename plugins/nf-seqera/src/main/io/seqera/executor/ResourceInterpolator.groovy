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

package io.seqera.executor

import java.util.regex.Pattern

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.util.MemoryUnit

/**
 * Rewrites a task command so that values derived from {@code task.memory} and
 * {@code task.cpus} become shell-variable placeholders, letting the scheduler fill
 * in the <em>actually allocated</em> value at launch time.
 *
 * <p>The problem (see sched#492): when resource optimization reduces a task's
 * allocation below its request, any value the command computed from the request
 * (e.g. a tool's {@code --max-mem} self-sizing flag) is now wrong and the task can
 * OOM. The command is rendered before the scheduler decides the allocation, so it
 * cannot be patched afterwards.
 *
 * <p>This works at materialization time, with no re-execution of the user script,
 * using three signals:
 * <ul>
 *   <li>the live command {@link GString} (already evaluated) — exact <em>location</em>
 *       of every interpolated value, recursing into nested GStrings;</li>
 *   <li>the raw script {@code source} — the interpolation <em>expressions</em>, aligned
 *       to the GString values by order;</li>
 *   <li>the AST variable references ({@code task} referenced at all) — a cheap gate.</li>
 * </ul>
 * A value is memory/cpus-derived when its source expression references
 * {@code task.memory}/{@code task.cpus} (directly, or through a local variable whose
 * definition does). Such numeric values are replaced with {@code ${SEQERA_TASK_MEM_n}} /
 * {@code ${SEQERA_TASK_CPUS_n}} and recorded as {@link Binding}s.
 *
 * <p>Safety: the parsed command template is validated against the live GString (the
 * literal text between interpolations must match, ignoring whitespace/escaping) before
 * any attribution, so a mis-identified template cannot produce wrong bindings. A
 * placeholder that would land inside shell single-quotes (where it would not expand) or
 * any parse/validation failure causes the original command to be returned unchanged —
 * never producing a worse command than today.
 *
 * @author Paolo Di Tommaso
 */
@Slf4j
@CompileStatic
class ResourceInterpolator {

    static enum Resource { MEMORY, CPUS }

    /** A placeholder to be resolved by the scheduler at launch time. */
    @Canonical
    static class Binding {
        String name
        Resource resource
        BigDecimal value
    }

    @Canonical
    static class Result {
        String script
        List<Binding> bindings
        boolean changed() { bindings != null && !bindings.isEmpty() }
    }

    private static class Template {
        List<String> exprs
        String literals
    }

    private static final String MEM_VAR = 'SEQERA_TASK_MEM_'
    private static final String CPUS_VAR = 'SEQERA_TASK_CPUS_'
    private static final Pattern PLACEHOLDER = Pattern.compile('\\$\\{SEQERA_TASK_')
    // anchored so task.memoryUsage / task.cpusPerNode / mytask.memory are NOT matched
    private static final Pattern MEM_REF = Pattern.compile(/(?<![A-Za-z0-9_.])task\s*\.\s*memory(?![A-Za-z0-9_])/)
    private static final Pattern CPUS_REF = Pattern.compile(/(?<![A-Za-z0-9_.])task\s*\.\s*cpus(?![A-Za-z0-9_])/)

    /**
     * @param command   the live command GString (retained by TaskRun before flattening)
     * @param source    the raw script-block source (TaskRun body source)
     * @param refNames  the AST-referenced variable names (TaskRun body valRefs)
     * @param reqCpus   the requested cpus the command was rendered with
     * @param reqMem    the requested memory the command was rendered with
     */
    static Result interpolate(GString command, String source, Set<String> refNames, int reqCpus, MemoryUnit reqMem) {
        final base = command != null ? command.toString() : null
        // -- gates: need a command GString, a source, a positive request, and a `task` reference
        if( command == null || !source || reqMem == null || reqCpus <= 0 )
            return new Result(base, [])
        if( refNames != null && !refNames.contains('task') )
            return new Result(base, [])

        try {
            // -- the command template = the last interpolating string literal in the source
            final templateText = lastStringLiteral(source)
            if( templateText == null )
                return new Result(base, [])
            final template = parseTemplate(templateText)
            // -- validate the parsed template actually corresponds to this command, so a
            //    mis-identified literal cannot produce wrong attributions
            if( !matches(template, command) )
                return new Result(base, [])

            final ctx = new Ctx(defs: parseGStringDefs(source))
            final placeholdered = render(command, template.exprs, ctx)

            // -- a placeholder that lands inside shell single-quotes would not expand: bail
            if( placeholderInSingleQuotes(placeholdered) ) {
                log.debug "[SEQERA] Resource interpolation skipped: placeholder inside single-quotes"
                return new Result(base, [])
            }
            return new Result(placeholdered, ctx.bindings)
        }
        catch( Exception e ) {
            log.debug "[SEQERA] Resource interpolation skipped: ${e.message}"
            return new Result(base, [])
        }
    }

    // ---- recursive render ---------------------------------------------------

    private static class Ctx {
        Map<String,String> defs
        List<Binding> bindings = []
        int memCount = 0
        int cpusCount = 0
    }

    /**
     * Walk a GString together with the ordered source expressions for its slots,
     * emitting resource-derived numeric values as shell-variable placeholders.
     */
    private static String render(GString g, List<String> exprs, Ctx ctx) {
        final strings = g.strings
        final values = g.values
        if( exprs == null || exprs.size() != values.length )
            throw new IllegalStateException("interpolation count mismatch")

        final sb = new StringBuilder()
        for( int i = 0; i < strings.length; i++ ) {
            sb.append(strings[i])
            if( i < values.length ) {
                final v = values[i]
                final e = exprs[i]?.trim()
                final kind = classify(e, ctx.defs)
                if( kind == Resource.MEMORY || kind == Resource.CPUS ) {
                    final num = asNumber(v)
                    if( num == null ) {
                        sb.append(str(v))
                    }
                    else {
                        final name = kind == Resource.MEMORY ? MEM_VAR + (ctx.memCount++) : CPUS_VAR + (ctx.cpusCount++)
                        ctx.bindings << new Binding(name, kind, num)
                        sb.append('${').append(name).append('}')
                    }
                }
                else if( kind == null && v instanceof GString && isLocalGString(e, ctx.defs) ) {
                    // one level of indirection: recurse into the local var's nested GString,
                    // but only if its definition template validates against the nested value
                    final nested = (GString) v
                    final defTpl = parseTemplate(ctx.defs[e])
                    if( matches(defTpl, nested) )
                        sb.append(render(nested, defTpl.exprs, ctx))
                    else
                        sb.append(str(v))
                }
                else {
                    sb.append(str(v))   // constant or joint expression -> unchanged
                }
            }
        }
        return sb.toString()
    }

    /**
     * @return MEMORY/CPUS when the expression references exactly one resource directly;
     *         null otherwise (constant, joint memory+cpus, or a local variable).
     */
    private static Resource classify(String expr, Map<String,String> defs) {
        if( !expr )
            return null
        final mem = MEM_REF.matcher(expr).find()
        final cpus = CPUS_REF.matcher(expr).find()
        if( mem && cpus )
            return null          // joint -> cannot scale by one resource (safe)
        if( mem )
            return Resource.MEMORY
        if( cpus )
            return Resource.CPUS
        return null
    }

    private static boolean isLocalGString(String expr, Map<String,String> defs) {
        if( !isLocalVar(expr) )
            return false
        final body = defs[expr]
        if( body == null )
            return false
        final mem = MEM_REF.matcher(body).find()
        final cpus = CPUS_REF.matcher(body).find()
        // recurse only for a single-resource local; joint locals are left untouched
        return mem ^ cpus
    }

    private static boolean isLocalVar(String expr) {
        expr != null && expr.matches(/[A-Za-z_][A-Za-z0-9_]*/)
    }

    // ---- template parsing & validation -------------------------------------

    /**
     * The parsed template literals must reproduce the live GString's literal segments,
     * and the interpolation count must match. Comparison ignores whitespace and
     * backslashes so source escaping / line-continuations do not cause false negatives
     * (which would only forgo the optimization, never mis-attribute).
     */
    private static boolean matches(Template t, GString g) {
        if( t == null || t.exprs.size() != g.values.length )
            return false
        return normalize(t.literals) == normalize(String.join('', g.strings))
    }

    private static String normalize(String s) {
        s == null ? '' : s.replaceAll(/[\s\\]/, '')
    }

    /**
     * Parse a GString template into its ordered interpolation expressions and the
     * concatenation of the literal text around them. Handles {@code ${ ... }} (balanced
     * braces) and {@code $dotted.path}; an escaped {@code \$} is treated as a literal.
     */
    static Template parseTemplate(String template) {
        final exprs = new ArrayList<String>()
        final lit = new StringBuilder()
        if( template != null ) {
            final n = template.length()
            int i = 0
            while( i < n ) {
                final c = template.charAt(i)
                if( c == ('\\' as char) && i + 1 < n ) {
                    lit.append(c).append(template.charAt(i + 1)); i += 2; continue
                }
                if( c == ('$' as char) && i + 1 < n ) {
                    final next = template.charAt(i + 1)
                    if( next == ('{' as char) ) {
                        int depth = 1, j = i + 2
                        while( j < n && depth > 0 ) {
                            final cj = template.charAt(j)
                            if( cj == ('{' as char) ) depth++
                            else if( cj == ('}' as char) ) { depth--; if( depth == 0 ) break }
                            j++
                        }
                        exprs << template.substring(i + 2, j)
                        i = j + 1
                        continue
                    }
                    else if( Character.isJavaIdentifierStart(next) ) {
                        int j = i + 1
                        while( j < n && (Character.isJavaIdentifierPart(template.charAt(j)) || template.charAt(j) == ('.' as char)) ) {
                            if( template.charAt(j) == ('.' as char) && (j + 1 >= n || !Character.isJavaIdentifierPart(template.charAt(j + 1))) )
                                break
                            j++
                        }
                        exprs << template.substring(i + 1, j)
                        i = j
                        continue
                    }
                }
                lit.append(c); i++
            }
        }
        final t = new Template()
        t.exprs = exprs
        t.literals = lit.toString()
        return t
    }

    /** @deprecated kept for tests; use {@link #parseTemplate}. */
    static List<String> parseInterpolations(String template) {
        parseTemplate(template).exprs
    }

    /**
     * Extract the inner text of the interpolating string literal that closes last in the
     * source — for a process script block this is the command template (the closure's
     * return value). Supports triple-double-quoted (""") and dollar-slashy ($/.../$);
     * other forms return null and are left unchanged (safe).
     */
    static String lastStringLiteral(String source) {
        String best = null
        int bestEnd = -1
        // triple-double-quoted
        final tdClose = source.lastIndexOf('"""')
        if( tdClose >= 0 ) {
            final open = source.lastIndexOf('"""', tdClose - 1)
            if( open >= 0 && open + 3 <= tdClose && tdClose + 3 > bestEnd ) {
                best = source.substring(open + 3, tdClose); bestEnd = tdClose + 3
            }
        }
        // dollar-slashy  $/ ... /$
        final dsClose = source.lastIndexOf('/$')
        if( dsClose >= 0 ) {
            final open = source.lastIndexOf('$/', dsClose - 1)
            if( open >= 0 && open + 2 <= dsClose && dsClose + 2 > bestEnd ) {
                best = source.substring(open + 2, dsClose); bestEnd = dsClose + 2
            }
        }
        return best
    }

    /**
     * Find {@code def NAME = <string literal>} definitions whose literal references a
     * resource, returning NAME -> the literal's inner template text. Used to follow one
     * level of local-variable indirection. The literal may sit inside a larger expression
     * (e.g. a ternary): the first resource-referencing literal is used.
     */
    static Map<String,String> parseGStringDefs(String source) {
        final defs = new LinkedHashMap<String,String>()
        final m = source =~ /(?m)\bdef\s+([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$/
        while( m.find() ) {
            final name = m.group(1) as String
            final rhs = m.group(2) as String
            final lit = firstResourceLiteral(rhs)
            if( lit != null )
                defs[name] = lit
        }
        return defs
    }

    private static String firstResourceLiteral(String rhs) {
        for( String q : ['"""', "'''", '"', "'"] ) {
            int from = 0
            while( true ) {
                final open = rhs.indexOf(q, from)
                if( open < 0 ) break
                final close = rhs.indexOf(q, open + q.length())
                if( close < 0 ) break
                final inner = rhs.substring(open + q.length(), close)
                if( MEM_REF.matcher(inner).find() || CPUS_REF.matcher(inner).find() )
                    return inner
                from = close + q.length()
            }
        }
        return null
    }

    // ---- shell safety -------------------------------------------------------

    /**
     * @return true if any {@code ${SEQERA_TASK_*}} placeholder falls inside a shell
     *         single-quoted region, where it would not be expanded at runtime.
     */
    static boolean placeholderInSingleQuotes(String s) {
        boolean inSingle = false, inDouble = false
        final n = s.length()
        for( int i = 0; i < n; i++ ) {
            final c = s.charAt(i)
            if( c == ('\\' as char) && !inSingle ) { i++; continue }   // backslash escape (not inside '...')
            if( c == ('\'' as char) && !inDouble ) { inSingle = !inSingle; continue }
            if( c == ('"' as char) && !inSingle ) { inDouble = !inDouble; continue }
            if( inSingle && s.startsWith('${SEQERA_TASK_', i) )
                return true
        }
        return false
    }

    // ---- helpers ------------------------------------------------------------

    private static String str(Object v) { v == null ? '' : v.toString() }

    private static BigDecimal asNumber(Object v) {
        if( v == null )
            return null
        try {
            return new BigDecimal(v.toString().trim())
        }
        catch( NumberFormatException e ) {
            return null
        }
    }
}
