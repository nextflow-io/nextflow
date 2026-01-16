/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.util

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * Helper methods to create class loaders with
 * additional classpaths.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ClassLoaderFactory {

    /**
     * Create a class loader with the given directories
     * added to the classpath.
     *
     * @param dirs
     */
    static GroovyClassLoader create(List<Path> dirs) {
        final gcl = new GroovyClassLoader()
        final libraries = resolveClassPaths(dirs)
        for( final lib : libraries ) {
            gcl.addClasspath(lib.complete().toString())
        }
        return gcl
    }

    /**
     * Given a list of directories, look for the files ending with
     * the '.jar' extension and return a list containing the original
     * directories and the JAR paths.
     *
     * @param dirs
     */
    static List<Path> resolveClassPaths(List<Path> dirs) {

        List<Path> result = []

        if( !dirs )
            return result

        for( final path : dirs ) {
            if( path.isFile() && path.name.endsWith('.jar') ) {
                result << path
            }
            else if( path.isDirectory() ) {
                result << path
                path.eachFileMatch( ~/.+\.jar$/ ) {
                    if( it.isFile() )
                        result << it
                }
            }
        }

        return result
    }

}
