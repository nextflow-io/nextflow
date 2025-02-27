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
 */

package nextflow.file

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import groovy.util.logging.Slf4j
import nextflow.plugin.Plugins
import org.pf4j.ExtensionPoint

/**
 * Define extension methods for supporting pluggable file remote file systems e.g. AWS S3 or Google Storage
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class FileSystemPathFactory implements ExtensionPoint {

    /**
     * Converts path uri string to the corresponding {@link Path} object
     *
     * @param uri
     *      A fully qualified path uri including the protocol scheme
     * @return
     *      A {@link Path} for the given path or {@code null} if protocol is unknown
     */
    abstract protected Path parseUri(String uri)

    /**
     * Converts a {@link Path} object to a fully qualified uri string
     *
     * @param path
     *      The {@link Path} object to be converted
     * @return
     *      The uri string corresponding to the specified path or {@code null}
     *      if the no provider is found
     */
    abstract protected String toUriString(Path path)

    /**
     * Bash helper library that implements the support for third-party storage
     * to be included in the job command wrapper script
     *
     * @return The Bash snippet implementing the support for third-party such as AWS S3
     * or {@code null} if not supported
     */
    abstract protected String getBashLib(Path target)

    /**
     * The name of a Bash helper function to upload a file to a remote file storage
     *
     * @return The name of the upload function or {@code null} if not supported
     */
    abstract protected String getUploadCmd(String source, Path target)


    static Path parse(String uri) {
        lookup0 { it.parseUri(uri) }
    }

    static String getUriString(Path path) {
        lookup0 { it.toUriString(path) }
    }

    static String bashLib(Path target) {
        lookup0 { it.getBashLib(target) }
    }

    static String uploadCmd(String source, Path target) {
        lookup0 { it.getUploadCmd(source, target) }
    }

    private static List<FileSystemPathFactory> factories0() {
        final factories = new ArrayList(10)
        final itr = Plugins.getPriorityExtensions(FileSystemPathFactory).iterator()
        while( itr.hasNext() )
            factories.add(itr.next())
        log.trace "File system path factories: ${factories}"
        return factories
    }

    private static <T> T lookup0( @ClosureParams(value = SimpleType, options = ['nextflow.file.FileSystemPathFactory']) Closure<T> criteria) {
        final factories = factories0()
        for( int i=0; i<factories.size(); i++ ) {
            final FileSystemPathFactory it = factories[i]
            final result = criteria.call(it)
            if( result!=null )
                return result
        }
        return null
    }

}
