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

package nextflow.file


import java.nio.file.Files
import java.nio.file.Path

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import org.yaml.snakeyaml.Yaml

/**
 * Add helper methods to {@link groovy.json.JsonSlurper} class
 *
 * Check resource file: `META-INF/services/org.codehaus.groovy.runtime.ExtensionModule
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class SlurperEx {

    /**
     * Parse a JSON data structure from content within a given File.
     *
     * @param file File containing JSON content
     * @return a data structure of lists and maps

     */
    static Object parse(JsonSlurper self, Path path) {
        return self.parse(Files.newInputStream(path))
    }

    /**
     * Parse a JSON data structure from content within a given File.
     *
     * @param file File containing JSON content
     * @param charset the charset for this File
     * @return a data structure of lists and maps
     */
    static Object parse(JsonSlurper self, Path path, String charset) {
        return self.parse(Files.newInputStream(path), charset)
    }

    /**
     * Load and parse a Yaml given the file {@link Path}
     *
     * @param self
     * @param path
     * @return
     */
    static Object load(Yaml self, Path path) {
        self.load(Files.newInputStream(path))
    }

    /**
     * Load and parse a Yaml given the file {@link File}
     *
     * @param self
     * @param file
     * @return
     */
    static Object load(Yaml self, File file) {
        self.load(Files.newInputStream(file.toPath()))
    }
}
