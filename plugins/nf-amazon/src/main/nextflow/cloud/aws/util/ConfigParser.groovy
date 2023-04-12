/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.cloud.aws.util

import java.nio.file.Path
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import org.apache.commons.lang.text.StrBuilder
/**
 * Parse and merge AWS config and credentials file
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ConfigParser {

    final private static Pattern KEY_VALUE = ~/\s*(\w+)\s*=.*/

    final Map<String, List<String>> content = new LinkedHashMap<>()

    ConfigParser parseConfig(Path path) {
        return parseConfig(path.text)
    }

    ConfigParser parseConfig(String text) {
        String current = null
        for( String line : text.readLines() ) {
            final section = parseSection(line)
            if( section ) {
                current = section
            }
            else if( current && line.trim() ) {
                final block = content.computeIfAbsent(current, (String it) -> new ArrayList<>())
                final key = findKey(line)
                final exists = key && block.any { findKey(it)==key }
                if( !key || !exists )
                    block.add(line)
            }
        }

        return this
    }

    protected String parseSection(String str) {
        def line = str.trim()
        if( !line.startsWith('[') || !line.endsWith(']') ) {
            return null
        }
        line = line.substring(1, line.size()-1)
        if( line.startsWith('profile '))
            line = line.substring('profile '.size())
        return line
    }

    String text() {
        final result = new StrBuilder()
        for( Map.Entry<String,List<String>> entry : content ) {
            result.append('[').append(entry.key).append(']\n')
            for( String line : entry.value ) {
                result.append(line).append('\n')
            }
        }
        return result.toString()
    }

    protected String findKey(String line) {
        final m = KEY_VALUE.matcher(line)
        return m.matches() ? m.group(1) : null
    }

}



