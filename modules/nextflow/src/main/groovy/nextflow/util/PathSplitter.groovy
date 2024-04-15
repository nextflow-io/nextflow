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

package nextflow.util

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import nextflow.file.FileHelper

/**
 * Split a path into two paths, the first component which may include the host name if it's a remote
 * url, and the remaining sub-components
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@EqualsAndHashCode
@CompileStatic
class PathSplitter {
    String head
    List<String> tail

    static PathSplitter parse(String path) {
        final baseUrl = FileHelper.baseUrl(path)
        if( !baseUrl )
            return split0(path, 0)

        if( baseUrl != path )
            return split0(path, baseUrl.length())
        else
            return new PathSplitter(baseUrl)
    }

    private static PathSplitter split0(String path, int offset) {
        def skip = path[offset]=='/' ? offset+1 : offset
        def start = path.indexOf('/', skip)
        if( start!=-1 && start<path.length() ) {
            def head = path.substring(0,start)
            def tail = path.substring(start)
            return new PathSplitter(head, tail.tokenize('/') ?: null)
        }
        else {
            return new PathSplitter(path)
        }
    }
}
