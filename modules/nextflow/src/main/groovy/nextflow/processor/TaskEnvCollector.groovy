/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.processor

import java.nio.file.Path

import groovy.transform.CompileStatic
/**
 * Implements the collection of environment variables
 * from the environment of a task execution.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class TaskEnvCollector {

    private Path workDir

    TaskEnvCollector(Path workDir) {
        this.workDir = workDir
    }

    Map collect() {
        final env = workDir.resolve(TaskRun.CMD_ENV).text
        final result = new HashMap(50)
        for( String line : env.readLines() ) {
            def tokens = tokenize0(line)
            def key = tokens[0]
            def value = tokens[1]
            if( !key )
                continue
            result.put(key, value)
        }
        return result
    }

    private List<String> tokenize0(String line) {
        int p = line.indexOf('=')
        return p == -1
                ? List.of(line, '')
                : List.of(line.substring(0,p), line.substring(p+1))
    }
}
