/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.splitter

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Defines file splitter extension methods
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SplitterEx {

    private static final Map EMPTY = Collections.emptyMap()

    static long countLines(Path self, Map opts=EMPTY) {
        new TextSplitter().options(opts) .target(self) .count()
    }

    static long countFasta(Path self, Map opts=EMPTY) {
        new FastaSplitter().options(opts) .target(self) .count()
    }

    static long countFastq(Path self, Map opts=EMPTY) {
        new FastqSplitter().options(opts) .target(self) .count()
    }

    static List splitText(Path self, Map opts=EMPTY) {
        new TextSplitter().options(opts) .target(self) .list()
    }

    static List splitFasta(Path self, Map opts=EMPTY) {
        new FastaSplitter().options(opts) .target(self) .list()
    }

    static List splitFastq(Path self, Map opts=EMPTY) {
        new FastqSplitter().options(opts) .target(self) .list()
    }

    static List splitCsv(Path self, Map opts=EMPTY) {
        new CsvSplitter().options(opts) .target(self) .list()
    }

}