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

import ch.qos.logback.classic.Level
import spock.lang.Specification

/**
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class LogFileFormatterTest extends Specification {

    private static final String ESC = ''

    // ----- parseEntryHeader -----

    def 'parseEntryHeader extracts level and capture positions from valid entry-start lines'() {
        given:
        final line = 'May-26 14:30:04.000 [main] ERROR nextflow.Session - Pipeline failed'

        when:
        final h = LogFileFormatter.parseEntryHeader(line)

        then:
        h != null
        h.level == Level.ERROR
        line.substring(h.timestampStart, h.timestampEnd) == 'May-26 14:30:04.000'
        line.substring(h.threadStart,    h.threadEnd)    == '[main]'
        line.substring(h.levelStart,     h.levelEnd)     == 'ERROR'
        line.substring(h.loggerStart,    h.loggerEnd)    == 'nextflow.Session'
        line.substring(h.dashStart,      h.dashEnd)      == '-'
        h.bodyStart > 0
    }

    def 'parseEntryHeader returns null for continuation lines'() {
        expect:
        LogFileFormatter.parseEntryHeader('\tat foo.Bar.baz(Bar.java:42)') == null
        LogFileFormatter.parseEntryHeader('java.lang.RuntimeException: oh no') == null
        LogFileFormatter.parseEntryHeader('') == null
        LogFileFormatter.parseEntryHeader('not a log line') == null
    }

    def 'parseEntryHeader detects every supported level'() {
        expect:
        LogFileFormatter.parseEntryHeader("May-26 14:30:00.000 [main] $word foo.Bar - msg").level == expected

        where:
        word    | expected
        'TRACE' | Level.TRACE
        'DEBUG' | Level.DEBUG
        'INFO'  | Level.INFO
        'WARN'  | Level.WARN
        'ERROR' | Level.ERROR
    }

    // ----- indicatorFor -----

    def 'indicatorFor returns plain 2 spaces when useColor is false'() {
        expect:
        LogFileFormatter.indicatorFor(level, false, true)  == '  '
        LogFileFormatter.indicatorFor(level, false, false) == '  '

        where:
        level << [Level.TRACE, Level.DEBUG, Level.INFO, Level.WARN, Level.ERROR, null]
    }

    def 'indicatorFor returns per-level block, with WARN/ERROR using 2 cells on entry-start and 1 cell on continuations'() {
        expect:
        LogFileFormatter.indicatorFor(Level.TRACE, true, true)  == LogFileFormatter.INDICATOR_TRACE
        LogFileFormatter.indicatorFor(Level.TRACE, true, false) == LogFileFormatter.INDICATOR_TRACE
        LogFileFormatter.indicatorFor(Level.DEBUG, true, true)  == LogFileFormatter.INDICATOR_PLAIN
        LogFileFormatter.indicatorFor(Level.DEBUG, true, false) == LogFileFormatter.INDICATOR_PLAIN
        LogFileFormatter.indicatorFor(Level.INFO,  true, true)  == LogFileFormatter.INDICATOR_INFO
        LogFileFormatter.indicatorFor(Level.INFO,  true, false) == LogFileFormatter.INDICATOR_INFO
        LogFileFormatter.indicatorFor(Level.WARN,  true, true)  == LogFileFormatter.INDICATOR_WARN_FILL
        LogFileFormatter.indicatorFor(Level.WARN,  true, false) == LogFileFormatter.INDICATOR_WARN
        LogFileFormatter.indicatorFor(Level.ERROR, true, true)  == LogFileFormatter.INDICATOR_ERROR_FILL
        LogFileFormatter.indicatorFor(Level.ERROR, true, false) == LogFileFormatter.INDICATOR_ERROR
        LogFileFormatter.indicatorFor(null,        true, true)  == LogFileFormatter.INDICATOR_PLAIN
    }

    def 'level indicator escape codes match the requested colors'() {
        expect:
        LogFileFormatter.INDICATOR_TRACE.startsWith("${ESC}[40m")        // black bg
        LogFileFormatter.INDICATOR_INFO  == LogFileFormatter.INDICATOR_TRACE   // INFO uses the same black cell
        LogFileFormatter.INDICATOR_WARN.startsWith("${ESC}[43m")         // yellow bg (continuation: 1 cell + plain space)
        LogFileFormatter.INDICATOR_WARN_FILL.startsWith("${ESC}[43m")    // yellow bg (entry-start: 2 cells, no trailing space)
        LogFileFormatter.INDICATOR_ERROR.startsWith("${ESC}[41m")        // red bg
        LogFileFormatter.INDICATOR_ERROR_FILL.startsWith("${ESC}[41m")
        LogFileFormatter.INDICATOR_PLAIN == '  '
        // entry-start fills end with reset (no plain trailing space) so bg flows continuously into timestamp
        LogFileFormatter.INDICATOR_WARN_FILL.endsWith("${ESC}[0m")
        LogFileFormatter.INDICATOR_ERROR_FILL.endsWith("${ESC}[0m")
    }

    // ----- format() -- no-color passes through -----

    def 'format returns line unchanged when useColor is false'() {
        given:
        final f = new LogFileFormatter(false)
        final line = 'May-26 14:30:00.000 [main] INFO  nextflow.Session - hello'
        final h = LogFileFormatter.parseEntryHeader(line)

        expect:
        f.format(line, Level.INFO, h) == line
    }

    // ----- format() -- TRACE/DEBUG wrapping -----

    def 'TRACE entry: indicator + dim wrapper, no body grammar applied'() {
        given:
        final f = new LogFileFormatter(true)
        final line = 'May-26 14:30:00.000 [main] TRACE nextflow.Session - msg'
        final h = LogFileFormatter.parseEntryHeader(line)

        when:
        final out = f.format(line, Level.TRACE, h)

        then:
        out.startsWith(LogFileFormatter.INDICATOR_TRACE)
        out.contains("${ESC}[2m")
        out.endsWith("${ESC}[0m")
        // body has no extra coloring (just the dim wrapper)
        !out.contains("${ESC}[36m")  // no cyan for thread
    }

    def 'DEBUG entry: plain prefix + dim wrapper AND grammar coloring'() {
        given:
        final f = new LogFileFormatter(true)
        final line = 'May-26 14:30:00.000 [main] DEBUG nextflow.Session - hello'
        final h = LogFileFormatter.parseEntryHeader(line)

        when:
        final out = f.format(line, Level.DEBUG, h)

        then:
        out.startsWith('  ')                  // DEBUG uses plain 2-space prefix
        out.startsWith(LogFileFormatter.INDICATOR_PLAIN + ESC + '[2m')
        out.endsWith("${ESC}[0m")
        // grammar styling still applied -- thread, logger colored -- and dim reactivated after each reset
        out.contains("${ESC}[36m[main]${ESC}[0m${ESC}[2m")     // cyan thread + reset + dim again
        out.contains("${ESC}[32mnextflow.Session${ESC}[0m${ESC}[2m")  // green logger + reset + dim again
    }

    def 'TRACE entry: dim wrapper, no grammar coloring'() {
        given:
        final f = new LogFileFormatter(true)
        final line = 'May-26 14:30:00.000 [main] TRACE nextflow.Session - msg'
        final h = LogFileFormatter.parseEntryHeader(line)

        when:
        final out = f.format(line, Level.TRACE, h)

        then:
        out == LogFileFormatter.INDICATOR_TRACE + "${ESC}[2m" + line + "${ESC}[0m"
        !out.contains("${ESC}[36m")  // no cyan thread coloring
    }

    def 'TRACE continuation lines inherit indicator and dim, no coloring'() {
        given:
        final f = new LogFileFormatter(true)
        final line = '\tat foo.Bar.baz(Bar.java:42)'

        when:
        final out = f.format(line, Level.TRACE, null)

        then:
        out == LogFileFormatter.INDICATOR_TRACE + "${ESC}[2m" + line + "${ESC}[0m"
    }

    def 'DEBUG continuation lines inherit indicator, dim AND grammar coloring'() {
        given:
        final f = new LogFileFormatter(true)
        final line = '\tat foo.Bar.baz(Bar.java:42)'

        when:
        final out = f.format(line, Level.DEBUG, null)

        then:
        out.startsWith('  ')
        out.contains("${ESC}[2m")
        // grammar still applies: 'at' keyword dim, frame name cyan
        out.contains("${ESC}[36mfoo.Bar.baz${ESC}[0m${ESC}[2m")
    }

    // ----- format() -- INFO/WARN/ERROR styling -----

    def 'INFO entry-start: single black-bg indicator + grammar styling, no level-token highlight'() {
        given:
        final f = new LogFileFormatter(true)
        final line = 'May-26 14:30:00.000 [main] INFO  nextflow.Session - hello'

        when:
        final out = f.format(line, Level.INFO, LogFileFormatter.parseEntryHeader(line))

        then:
        out.startsWith(LogFileFormatter.INDICATOR_INFO)
        !out.contains("${ESC}[1mINFO${ESC}[0m")
        out.contains(' INFO  ')
        // grammar applies normally to INFO: dim timestamp, cyan thread, green logger
        out.contains("${ESC}[2mMay-26 14:30:00.000${ESC}[0m")
        out.contains("${ESC}[36m[main]${ESC}[0m")
        out.contains("${ESC}[32mnextflow.Session${ESC}[0m")
    }

    def 'WARN entry-start: 2-cell indicator + continuous yellow bg through level token'() {
        given:
        final f = new LogFileFormatter(true)
        final line = 'May-26 14:30:00.000 [main] WARN  nextflow.Session - heads up'

        when:
        final out = f.format(line, Level.WARN, LogFileFormatter.parseEntryHeader(line))

        then:
        out.startsWith(LogFileFormatter.INDICATOR_WARN_FILL)
        out.contains("${ESC}[30;43mMay-26 14:30:00.000 [main] WARN${ESC}[0m")
        !out.contains("${ESC}[36m[main]${ESC}[0m")
        out.contains("${ESC}[0m${ESC}[33m")
        out.contains("${ESC}[32mnextflow.Session${ESC}[0m")
    }

    def 'ERROR entry-start: 2-cell indicator + continuous red bg through level token'() {
        given:
        final f = new LogFileFormatter(true)
        final line = 'May-26 14:30:00.000 [main] ERROR nextflow.Session - oh no'

        when:
        final out = f.format(line, Level.ERROR, LogFileFormatter.parseEntryHeader(line))

        then:
        out.startsWith(LogFileFormatter.INDICATOR_ERROR_FILL)
        out.contains("${ESC}[37;41mMay-26 14:30:00.000 [main] ERROR${ESC}[0m")
        !out.contains("${ESC}[36m[main]${ESC}[0m")
        out.contains("${ESC}[0m${ESC}[31m")
        out.contains("${ESC}[32mnextflow.Session${ESC}[0m")
    }

    def 'WARN continuation line: single yellow cell + yellow body fg'() {
        given:
        final f = new LogFileFormatter(true)
        final line = '   some extra warning context'

        when:
        final out = f.format(line, Level.WARN, null)

        then:
        out.startsWith(LogFileFormatter.INDICATOR_WARN + "${ESC}[33m")
        !out.startsWith(LogFileFormatter.INDICATOR_WARN_FILL)   // continuations use 1 cell, not 2
        out.endsWith("${ESC}[0m")
    }

    def 'ERROR continuation line: single red cell + red body fg + grammar coloring still applies'() {
        given:
        final f = new LogFileFormatter(true)
        final line = '\tat foo.Bar.baz(Bar.java:42)'

        when:
        final out = f.format(line, Level.ERROR, null)

        then:
        out.startsWith(LogFileFormatter.INDICATOR_ERROR + "${ESC}[31m")
        !out.startsWith(LogFileFormatter.INDICATOR_ERROR_FILL)
        out.contains("${ESC}[36mfoo.Bar.baz${ESC}[0m")
        out.contains("${ESC}[0m${ESC}[31m")
    }

    // ----- format() -- grammar body rules -----

    def 'format highlights work-hash in entry-start body'() {
        given:
        final f = new LogFileFormatter(true)
        final line = 'May-26 14:30:00.000 [main] INFO  nextflow.Session - [12/abcdef] Submitting task'
        final h = LogFileFormatter.parseEntryHeader(line)

        when:
        final out = f.format(line, Level.INFO, h)

        then:
        out.contains("${ESC}[33m[12/abcdef]${ESC}[0m")
    }

    def 'format highlights stack-trace continuation lines'() {
        given:
        final f = new LogFileFormatter(true)
        final continuation = '\tat foo.Bar.baz(Bar.java:42)'

        when:
        final out = f.format(continuation, Level.ERROR, null)

        then:
        out.startsWith(LogFileFormatter.INDICATOR_ERROR)              // 1-cell continuation indicator
        !out.startsWith(LogFileFormatter.INDICATOR_ERROR_FILL)
        out.contains("${ESC}[2mat${ESC}[0m")
        out.contains("${ESC}[36mfoo.Bar.baz${ESC}[0m")
    }

    def 'format highlights exception lines'() {
        given:
        final f = new LogFileFormatter(true)
        final line = 'java.lang.RuntimeException: oh no'

        when:
        final out = f.format(line, Level.ERROR, null)

        then:
        out.contains("${ESC}[1;31mjava.lang.RuntimeException${ESC}[0m")
        out.contains("${ESC}[31moh no${ESC}[0m")
    }

    def 'format highlights URLs as underlined blue'() {
        given:
        final f = new LogFileFormatter(true)
        final line = 'May-26 14:30:00.000 [main] INFO  nextflow.Session - see https://nextflow.io/docs for help'
        final h = LogFileFormatter.parseEntryHeader(line)

        when:
        final out = f.format(line, Level.INFO, h)

        then:
        out.contains("${ESC}[4;34mhttps://nextflow.io/docs${ESC}[0m")
    }

    def 'format leaves pre-header lines (no level) unstyled apart from plain prefix'() {
        given:
        final f = new LogFileFormatter(true)
        final line = 'some preamble before any log entry'

        when:
        final out = f.format(line, null, null)

        then:
        out == '  ' + line  // plain 2-space prefix, no styling
    }

    // ----- format() -- useColor=false path adds no prefix -----

    def 'useColor=false preserves bytes -- no prefix, no ANSI'() {
        given:
        final f = new LogFileFormatter(false)
        final line = 'May-26 14:30:00.000 [main] INFO  nextflow.Session - hi'

        expect:
        f.format(line, Level.INFO, LogFileFormatter.parseEntryHeader(line)) == line
        f.format(line, Level.TRACE, LogFileFormatter.parseEntryHeader(line)) == line
        f.format('continuation', Level.ERROR, null) == 'continuation'
    }
}
