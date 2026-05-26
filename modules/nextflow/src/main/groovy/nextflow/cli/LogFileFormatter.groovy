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

import java.util.regex.Matcher
import java.util.regex.Pattern

import ch.qos.logback.classic.Level
import groovy.transform.CompileStatic
import groovy.transform.PackageScope

/**
 * Renders Nextflow log file lines for terminal display.
 *
 * Each line is prefixed with a 2-character level indicator and the message body is
 * highlighted using rules ported from the nextflow-log TextMate grammar shipped with
 * vscode-language-nextflow.
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@CompileStatic
class LogFileFormatter {

    static final String ESC = ''
    static final String ANSI_RESET = ESC + '[0m'
    static final String ANSI_DIM   = ESC + '[2m'

    // Level indicator: 2-char prefix. WARN/ERROR continuations use 1 colored cell + plain
    // space; entry-start lines use 2 colored cells so the bg flows continuously into the
    // timestamp emphasis with no visible gap.
    static final String INDICATOR_TRACE      = ESC + '[40m ' + ANSI_RESET + ' '
    static final String INDICATOR_INFO       = INDICATOR_TRACE
    static final String INDICATOR_WARN       = ESC + '[43m ' + ANSI_RESET + ' '
    static final String INDICATOR_WARN_FILL  = ESC + '[43m  ' + ANSI_RESET
    static final String INDICATOR_ERROR      = ESC + '[41m ' + ANSI_RESET + ' '
    static final String INDICATOR_ERROR_FILL = ESC + '[41m  ' + ANSI_RESET
    static final String INDICATOR_PLAIN      = '  '

    static final String WARN_HIGHLIGHT  = ESC + '[30;43m'   // black fg on yellow bg
    static final String ERROR_HIGHLIGHT = ESC + '[37;41m'   // white fg on red bg
    static final String WARN_BODY_FG    = ESC + '[33m'      // yellow fg
    static final String ERROR_BODY_FG   = ESC + '[31m'      // red fg

    // Scopes referenced by stylizeHeader (separate from the body grammar rules)
    private static final String SCOPE_TIMESTAMP = 'constant.other.timestamp.nextflow-log'
    private static final String SCOPE_THREAD    = 'entity.name.tag.thread.nextflow-log'
    private static final String SCOPE_LOGGER    = 'entity.name.class.logger.nextflow-log'
    private static final String SCOPE_SEPARATOR = 'punctuation.separator.nextflow-log'

    /** Entry-start pattern. Captures: 1=timestamp, 2=thread, 3=level, 4=logger, 5=dash. */
    static final Pattern ENTRY_HEADER = ~/^([A-Z][a-z]{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3}) (\[[^\]]+\]) (TRACE|DEBUG|INFO|WARN|ERROR) +([\w$.]+) (-)/

    /** Captured offsets of an entry-start match. Immutable snapshot — does not keep the Matcher. */
    @PackageScope
    @CompileStatic
    static class EntryHeader {
        final Level level
        final int timestampStart, timestampEnd
        final int threadStart,    threadEnd
        final int levelStart,     levelEnd
        final int loggerStart,    loggerEnd
        final int dashStart,      dashEnd
        final int bodyStart

        EntryHeader(Matcher m) {
            this.level = Level.toLevel(m.group(3))
            this.timestampStart = m.start(1); this.timestampEnd = m.end(1)
            this.threadStart    = m.start(2); this.threadEnd    = m.end(2)
            this.levelStart     = m.start(3); this.levelEnd     = m.end(3)
            this.loggerStart    = m.start(4); this.loggerEnd    = m.end(4)
            this.dashStart      = m.start(5); this.dashEnd      = m.end(5)
            this.bodyStart      = m.end()
        }
    }

    static EntryHeader parseEntryHeader(String line) {
        final m = ENTRY_HEADER.matcher(line)
        m.find() ? new EntryHeader(m) : null
    }

    /** Mapping from TextMate-style scope name to ANSI escape sequence. */
    static final Map<String, String> THEME = Collections.unmodifiableMap([
        (SCOPE_TIMESTAMP)                               : ESC + '[2m',
        (SCOPE_THREAD)                                  : ESC + '[36m',
        (SCOPE_LOGGER)                                  : ESC + '[32m',
        (SCOPE_SEPARATOR)                               : ESC + '[2m',
        'keyword.control.section.nextflow-log'          : ESC + '[1;33m',
        'entity.name.type.exception.nextflow-log'       : ESC + '[1;31m',
        'string.unquoted.exception-message.nextflow-log': ESC + '[31m',
        'keyword.control.at.nextflow-log'               : ESC + '[2m',
        'entity.name.function.stack-frame.nextflow-log' : ESC + '[36m',
        'punctuation.section.parens.begin.nextflow-log' : ESC + '[2m',
        'punctuation.section.parens.end.nextflow-log'   : ESC + '[2m',
        'constant.other.source-location.nextflow-log'   : ESC + '[2m',
        'constant.numeric.hash.work.nextflow-log'       : ESC + '[33m',
        'string.unquoted.path.work.nextflow-log'        : ESC + '[36m',
        'keyword.other.process-marker.nextflow-log'     : ESC + '[35m',
        'entity.name.function.process.nextflow-log'     : ESC + '[1;36m',
        'meta.process-instance.nextflow-log'            : ESC + '[2m',
        'variable.parameter.port.nextflow-log'          : ESC + '[34m',
        'support.type.channel-type.nextflow-log'        : ESC + '[35m',
        'constant.language.state.nextflow-log'          : ESC + '[1;32m',
        'variable.parameter.nextflow-log'               : ESC + '[34m',
        'variable.other.channel.nextflow-log'           : ESC + '[36m',
        'punctuation.separator.key-value.nextflow-log'  : ESC + '[2m',
        'support.constant.script-id.nextflow-log'       : ESC + '[33m',
        'support.module.plugin.nextflow-log'            : ESC + '[35m',
        'punctuation.separator.version.nextflow-log'    : ESC + '[2m',
        'support.constant.version.nextflow-log'         : ESC + '[36m',
        'markup.underline.link.nextflow-log'            : ESC + '[4;34m',
        'string.quoted.double.nextflow-log'             : ESC + '[32m',
        'string.quoted.single.nextflow-log'             : ESC + '[32m',
        'string.unquoted.path.nextflow-log'             : ESC + '[36m',
        'string.unquoted.value.nextflow-log'            : ESC + '[32m',
        'constant.numeric.nextflow-log'                 : ESC + '[33m',
        'keyword.other.unit.nextflow-log'               : ESC + '[2;33m',
    ] as Map<String, String>)

    /**
     * A body grammar rule. Capture group indexes are sorted ascending and their ANSI codes
     * resolved from {@link #THEME} at construction so apply-time is allocation-free.
     */
    @PackageScope
    @CompileStatic
    static class BodyRule {
        final Pattern pattern
        final int[] groupIndexes
        final String[] groupCodes

        BodyRule(Pattern pattern, Map<Integer, String> captures) {
            this.pattern = pattern
            final ids = new ArrayList<Integer>(captures.keySet())
            Collections.sort(ids)
            groupIndexes = new int[ids.size()]
            groupCodes   = new String[ids.size()]
            for( int i = 0; i < ids.size(); i++ ) {
                final id = ids.get(i)
                groupIndexes[i] = id
                groupCodes[i]   = THEME.get(captures.get(id))
            }
        }
    }

    /** Body patterns, in priority order — `message-body` rules from the grammar. */
    static final List<BodyRule> BODY_RULES = Collections.unmodifiableList([
        new BodyRule(
            ~/^(Caused by|Command executed|Command exit status|Command output|Command error|Work dir|Tip):/,
            [(1): 'keyword.control.section.nextflow-log']),
        new BodyRule(
            ~/^([\w.$]+(?:Exception|Error)): (.+)$/,
            [(1): 'entity.name.type.exception.nextflow-log',
             (2): 'string.unquoted.exception-message.nextflow-log']),
        new BodyRule(
            ~/^(\t+)(at) ([\w$.\/<>]+)(\()([^)]+)(\))/,
            [(2): 'keyword.control.at.nextflow-log',
             (3): 'entity.name.function.stack-frame.nextflow-log',
             (4): 'punctuation.section.parens.begin.nextflow-log',
             (5): 'constant.other.source-location.nextflow-log',
             (6): 'punctuation.section.parens.end.nextflow-log']),
        new BodyRule(
            ~/\[[a-f0-9]{2}\/[a-f0-9]{6,}\]/,
            [(0): 'constant.numeric.hash.work.nextflow-log']),
        new BodyRule(
            ~/\/work\/[a-f0-9]{2}\/[a-f0-9]{6,}[\w.\/-]*/,
            [(0): 'string.unquoted.path.work.nextflow-log']),
        new BodyRule(
            ~/(process|workflow)( > )([A-Za-z0-9_:]+)((?:\s*\([^)]+\))?)/,
            [(1): 'keyword.other.process-marker.nextflow-log',
             (2): 'punctuation.separator.nextflow-log',
             (3): 'entity.name.function.process.nextflow-log',
             (4): 'meta.process-instance.nextflow-log']),
        new BodyRule(
            ~/^(\[process\]) ([A-Za-z0-9_:]+)/,
            [(1): 'keyword.other.process-marker.nextflow-log',
             (2): 'entity.name.function.process.nextflow-log']),
        new BodyRule(
            ~/(port \d+): (\((value|queue|cntrl)\)) (\S+)\s*; (channel): (.+)$/,
            [(1): 'variable.parameter.port.nextflow-log',
             (2): 'support.type.channel-type.nextflow-log',
             (4): 'constant.language.state.nextflow-log',
             (5): 'variable.parameter.nextflow-log',
             (6): 'variable.other.channel.nextflow-log']),
        new BodyRule(
            ~/^( +)(status)(=)(\S+)/,
            [(2): 'variable.parameter.nextflow-log',
             (3): 'punctuation.separator.key-value.nextflow-log',
             (4): 'constant.language.state.nextflow-log']),
        new BodyRule(
            ~/\bScript_[a-f0-9]+\b/,
            [(0): 'support.constant.script-id.nextflow-log']),
        new BodyRule(
            ~/\b([a-zA-Z][\w-]*)(@)(\d[\d.]*)\b/,
            [(1): 'support.module.plugin.nextflow-log',
             (2): 'punctuation.separator.version.nextflow-log',
             (3): 'support.constant.version.nextflow-log']),
        new BodyRule(
            ~/https?:\/\/[^\s)]+/,
            [(0): 'markup.underline.link.nextflow-log']),
        new BodyRule(
            ~/"[^"\n]*"/,
            [(0): 'string.quoted.double.nextflow-log']),
        new BodyRule(
            ~/'[^'\n]*'/,
            [(0): 'string.quoted.single.nextflow-log']),
        new BodyRule(
            ~/(?<![\w.\/@-])\/(?:[\w.-]+\/)*[\w.-]+\/?/,
            [(0): 'string.unquoted.path.nextflow-log']),
        new BodyRule(
            ~/\b([a-zA-Z_][\w.-]*)(=)([^\s;,]+)/,
            [(1): 'variable.parameter.nextflow-log',
             (2): 'punctuation.separator.key-value.nextflow-log',
             (3): 'string.unquoted.value.nextflow-log']),
        new BodyRule(
            ~/\b(\d+(?:\.\d+)?)(\s*(?:ms|s|m|h|d|GB|MB|KB|TB|B|%))?\b/,
            [(1): 'constant.numeric.nextflow-log',
             (2): 'keyword.other.unit.nextflow-log']),
    ] as List<BodyRule>)

    final boolean useColor

    LogFileFormatter(boolean useColor) {
        this.useColor = useColor
    }

    static String indicatorFor(Level level, boolean useColor, boolean isEntryStart) {
        if( !useColor || level == null )
            return INDICATOR_PLAIN
        if( level == Level.TRACE ) return INDICATOR_TRACE
        if( level == Level.INFO )  return INDICATOR_INFO
        if( level == Level.WARN )  return isEntryStart ? INDICATOR_WARN_FILL  : INDICATOR_WARN
        if( level == Level.ERROR ) return isEntryStart ? INDICATOR_ERROR_FILL : INDICATOR_ERROR
        return INDICATOR_PLAIN  // DEBUG and anything else
    }

    /**
     * Render a log line for output. {@code header} is non-null for entry-start lines and
     * null for continuation lines (which inherit {@code currentLevel} from their entry).
     */
    String format(String line, Level currentLevel, EntryHeader header) {
        if( !useColor )
            return line
        final indicator = indicatorFor(currentLevel, true, header != null)
        if( currentLevel == Level.TRACE )
            return indicator + ANSI_DIM + line + ANSI_RESET
        if( currentLevel == Level.DEBUG )
            return indicator + persistStyle(stylize(line, header, null), ANSI_DIM)
        if( currentLevel == Level.WARN )
            return indicator + persistStyle(stylize(line, header, WARN_HIGHLIGHT), WARN_BODY_FG)
        if( currentLevel == Level.ERROR )
            return indicator + persistStyle(stylize(line, header, ERROR_HIGHLIGHT), ERROR_BODY_FG)
        return indicator + stylize(line, header, null)
    }

    /**
     * Wrap a styled span in a persistent attribute. Each internal {@code [0m} reset is
     * followed by a fresh copy of {@code style} so the attribute remains active across the
     * colored tokens produced by the grammar.
     */
    private String persistStyle(String styled, String style) {
        style + styled.replace(ANSI_RESET, ANSI_RESET + style) + ANSI_RESET
    }

    private String stylize(String line, EntryHeader header, String headerEmphasis) {
        if( header == null )
            return tokenizeBody(line)
        return stylizeHeader(line, header, headerEmphasis) + tokenizeBody(line.substring(header.bodyStart))
    }

    private String stylizeHeader(String line, EntryHeader h, String headerEmphasis) {
        final sb = new StringBuilder()
        if( headerEmphasis != null ) {
            // WARN/ERROR: emphasize the whole header region (timestamp + thread + level) as
            // a single solid bg block so the indicator's bg flows continuously across.
            sb.append(headerEmphasis).append(line, 0, h.levelEnd).append(ANSI_RESET)
        }
        else {
            sb.append(line, 0, h.timestampStart)
            appendColored(sb, line, h.timestampStart, h.timestampEnd, SCOPE_TIMESTAMP)
            sb.append(line, h.timestampEnd, h.threadStart)
            appendColored(sb, line, h.threadStart, h.threadEnd, SCOPE_THREAD)
            sb.append(line, h.threadEnd, h.levelStart)
            // level token plain — indicator block carries the signal
            sb.append(line, h.levelStart, h.levelEnd)
        }
        sb.append(line, h.levelEnd, h.loggerStart)
        appendColored(sb, line, h.loggerStart, h.loggerEnd, SCOPE_LOGGER)
        sb.append(line, h.loggerEnd, h.dashStart)
        appendColored(sb, line, h.dashStart, h.dashEnd, SCOPE_SEPARATOR)
        return sb.toString()
    }

    private static void appendColored(StringBuilder sb, CharSequence src, int start, int end, String scope) {
        final code = THEME.get(scope)
        if( code )
            sb.append(code).append(src, start, end).append(ANSI_RESET)
        else
            sb.append(src, start, end)
    }

    private String tokenizeBody(String body) {
        if( body.isEmpty() )
            return body
        final n = body.length()
        final sb = new StringBuilder()
        int pos = 0
        while( pos < n ) {
            BodyRule bestRule = null
            Matcher bestMatcher = null
            int bestStart = Integer.MAX_VALUE
            for( BodyRule rule : BODY_RULES ) {
                final m = rule.pattern.matcher(body)
                if( m.find(pos) && m.start() < bestStart ) {
                    bestStart = m.start()
                    bestRule = rule
                    bestMatcher = m
                }
            }
            if( bestRule == null ) {
                sb.append(body, pos, n)
                break
            }
            if( pos < bestStart )
                sb.append(body, pos, bestStart)
            applyRule(sb, body, bestMatcher, bestRule)
            // guard against zero-length matches that would loop forever
            pos = bestMatcher.end() > bestStart ? bestMatcher.end() : bestStart + 1
        }
        return sb.toString()
    }

    private static void applyRule(StringBuilder sb, String body, Matcher m, BodyRule rule) {
        final indexes = rule.groupIndexes
        final codes   = rule.groupCodes
        // whole-match coloring
        if( indexes.length == 1 && indexes[0] == 0 ) {
            final code = codes[0]
            if( code )
                sb.append(code).append(body, m.start(), m.end()).append(ANSI_RESET)
            else
                sb.append(body, m.start(), m.end())
            return
        }
        // per-group coloring; preserve unmarked text between captured groups
        int cursor = m.start()
        for( int i = 0; i < indexes.length; i++ ) {
            final g = indexes[i]
            final s = m.start(g)
            final e = m.end(g)
            if( s < 0 )
                continue
            if( s > cursor )
                sb.append(body, cursor, s)
            final code = codes[i]
            if( code )
                sb.append(code).append(body, s, e).append(ANSI_RESET)
            else
                sb.append(body, s, e)
            cursor = e
        }
        if( cursor < m.end() )
            sb.append(body, cursor, m.end())
    }
}
