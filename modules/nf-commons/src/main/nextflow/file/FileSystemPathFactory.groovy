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

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
/**
 * Generic interface
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
abstract class FileSystemPathFactory {

    abstract protected Path parseUri(String uri)

    static Path parse(String uri) {
        def factories = factories0()
        for( int i=0; i<factories.size(); i++ ) {
            final result = factories[i].parseUri(uri)
            if( result )
                return result
        }
        return null
    }

    @Memoized
    private static List<FileSystemPathFactory> factories0() {
        final result = new ArrayList(10)
        final itr = ServiceLoader.load(FileSystemPathFactory).iterator()
        while( itr.hasNext() )
            result.add(itr.next())
        return result
    }
}
