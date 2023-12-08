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

package nextflow.script

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.IllegalFileException
import nextflow.file.FilePatternSplitter
import nextflow.util.BlankSeparatedList
/**
 * Models a process file output, which defines a file
 * or set of files to be unstaged from a task work directory.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ProcessFileOutput implements PathArityAware {

    private Object target

    /**
     * When true it will not fail if no files are found.
     */
    boolean optional

    /**
     * When true it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: true)
     */
    boolean followLinks = true

    /**
     * When true the specified name is interpreted as a glob pattern (default: true)
     */
    boolean glob = true

    /**
     * When {@code true} star wildcard (*) matches hidden files (files starting with a dot char)
     * By default it does not, coherently with linux bash rule
     */
    boolean hidden

    /**
     * When {@code true} file pattern includes input files as well as output files.
     * By default a file pattern matches only against files produced by the process, not
     * the ones received as input
     */
    boolean includeInputs

    /**
     * Maximum number of directory levels to visit (default: no limit)
     */
    Integer maxDepth

    /**
     * The type of path to output, either 'file', 'dir' or 'any'
     */
    String type

    ProcessFileOutput(Object target, Map opts) {
        this.target = target

        for( Map.Entry<String,?> entry : opts )
            setProperty(entry.key, entry.value)
    }

    List<String> getFilePatterns(Map context, Path workDir) {
        final entry = resolve(context, target)

        if( !entry )
            return []

        // -- single path
        if( entry instanceof Path )
            return [ relativize(entry, workDir) ]

        // -- multiple paths
        if( entry instanceof BlankSeparatedList || entry instanceof List )
            return entry.collect( path -> relativize(path.toString(), workDir) )

        // -- literal file name
        return [ relativize(entry.toString(), workDir) ]
    }

    protected Object resolve(Map ctx, Object value) {

        if( value instanceof GString )
            return value.cloneAsLazy(ctx)

        if( value instanceof Closure )
            return ctx.with(value)

        return value.toString()
    }

    protected String relativize(String path, Path workDir) {
        if( !path.startsWith('/') )
            return path

        final dir = workDir.toString()
        if( !path.startsWith(dir) )
            throw new IllegalFileException("File `$path` is outside the scope of the process work directory: $workDir")

        if( path.length()-dir.length()<2 )
            throw new IllegalFileException("Missing output file name")

        return path.substring(dir.size()+1)
    }

    protected String relativize(Path path, Path workDir) {
        if( !path.isAbsolute() )
            return glob ? FilePatternSplitter.GLOB.escape(path) : path

        if( !path.startsWith(workDir) )
            throw new IllegalFileException("File `$path` is outside the scope of the process work directory: $workDir")

        if( path.nameCount == workDir.nameCount )
            throw new IllegalFileException("Missing output file name")

        final rel = path.subpath(workDir.getNameCount(), path.getNameCount())
        return glob ? FilePatternSplitter.GLOB.escape(rel) : rel
    }
}
