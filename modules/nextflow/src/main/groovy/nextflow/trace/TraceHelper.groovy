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
 *
 */

package nextflow.trace

import java.nio.charset.Charset
import java.nio.file.FileAlreadyExistsException
import java.nio.file.Files
import java.nio.file.OpenOption
import java.nio.file.Path
import static java.nio.file.StandardOpenOption.*
import java.time.Instant
import java.time.ZoneId
import java.time.format.DateTimeFormatter

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.exception.AbortOperationException

/**
 * Helper methods for trace and observer feature
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TraceHelper {

    /* only for testing -- do not use */
    static protected String testTimestampFmt

    static public final Instant launchTimestamp = Instant.now()

    static String launchTimestampFmt() {
        return testTimestampFmt ?: fmt(launchTimestamp)
    }

    @Memoized
    static private String fmt(Instant instance) {
        final formatter = DateTimeFormatter.ofPattern("yyyyMMdd-A").withZone(ZoneId.systemDefault())
        return formatter.format(instance)
    }

    static private OpenOption[] openOptions(Boolean overwrite) {
        overwrite ? [ CREATE, WRITE, TRUNCATE_EXISTING ] as OpenOption[] : [ CREATE_NEW, WRITE ] as OpenOption[]
    }

    static BufferedWriter newFileWriter(Path path, boolean overwrite, String type) {
        try {
            Files.newBufferedWriter(path, Charset.defaultCharset(), openOptions(overwrite))
        }
        catch (FileAlreadyExistsException e) {
            throw new AbortOperationException("$type file already exists: ${path.toUriString()} -- enable the '${type.toLowerCase()}.overwrite' option in your config file to overwrite existing files", e)
        }
    }
}
