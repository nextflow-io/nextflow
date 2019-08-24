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

package nextflow.script.params

import java.nio.file.Path

import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.IllegalFileException
import nextflow.file.FilePatternSplitter
import nextflow.script.TokenVar
import nextflow.util.BlankSeparatedList
/**
 * Model a process *file* output parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
class FileOutParam extends BaseOutParam implements OutParam, OptionalParam, PathQualifier {

    /**
     * ONLY FOR TESTING DO NOT USE
     */
    protected FileOutParam(Map params) {
        super(new Binding(), [])
    }

    /**
     * The character used to separate multiple names (pattern) in the output specification
     *
     * This is only used by `file` qualifier. It's not supposed to be used anymore
     * by the new `path` qualifier.
     *
     */
    @Deprecated
    String separatorChar = ':'

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
     * The type of path to output, either {@code file}, {@code dir} or {@code any}
     */
    String type

    /**
     * Maximum number of directory levels to visit (default: no limit)
     */
    Integer maxDepth

    /**
     * When true it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: true)
     */
    boolean followLinks = true

    boolean glob = true

    private GString gstring
    private Closure<String> dynamicObj
    private String filePattern
    private boolean pathQualifier

    /**
     * @return {@code true} when the file name is parametric i.e contains a variable name to be resolved, {@code false} otherwise
     */
    boolean isDynamic() { dynamicObj || gstring != null }

    @Deprecated
    FileOutParam separatorChar( String value ) {
        this.separatorChar = value
        return this
    }

    @Deprecated
    FileOutParam includeInputs( boolean flag ) {
        this.includeInputs = flag
        return this
    }

    @Deprecated
    FileOutParam includeHidden( boolean flag ) {
        this.hidden = flag
        return this
    }

    @Deprecated
    FileOutParam hidden( boolean flag ) {
        this.hidden = flag
        return this
    }

    @Deprecated
    FileOutParam type( String value ) {
        assert value in ['file','dir','any']
        type = value
        return this
    }

    @Deprecated
    FileOutParam maxDepth( int value ) {
        maxDepth = value
        return this
    }

    @Deprecated
    FileOutParam followLinks( boolean value ) {
        followLinks = value
        return this
    }

    @Deprecated
    FileOutParam glob( boolean value ) {
        glob = value
        return this
    }

    @Override
    BaseOutParam bind( obj ) {

        if( obj instanceof GString ) {
            gstring = obj
            return this
        }

        if( obj instanceof TokenVar ) {
            this.nameObj = obj.name
            dynamicObj = { delegate.containsKey(obj.name) ? delegate.get(obj.name): obj.name }
            return this
        }

        if( obj instanceof Closure ) {
            dynamicObj = obj
            return this
        }

        this.filePattern = obj.toString()
        return this
    }

    List<String> getFilePatterns(Map context, Path workDir) {

        def entry = null
        if( dynamicObj ) {
            entry = context.with(dynamicObj)
        }
        else if( gstring != null ) {
            def strict = (getName() == null)
            try {
                entry = gstring.cloneAsLazy(context)
            }
            catch( MissingPropertyException e ) {
                if( strict )
                    throw e
            }
        }
        else {
            entry = filePattern
        }

        if( !entry )
            return []

        if( entry instanceof Path )
            return [ relativize(entry, workDir) ]

        // handle a collection of files
        if( entry instanceof BlankSeparatedList || entry instanceof List ) {
            return entry.collect { relativize(it.toString(), workDir) }
        }

        // normalize to a string object
        final nameString = entry.toString()
        if( separatorChar && nameString.contains(separatorChar) ) {
            return nameString.split(/\${separatorChar}/).collect { String it-> relativize(it, workDir) }
        }

        return [relativize(nameString, workDir)]

    }

    @PackageScope String getFilePattern() { filePattern }

    @PackageScope
    static String clean(String path) {
        while (path.startsWith('/') ) {
            path = path.substring(1)
        }
        return path
    }

    @PackageScope
    String relativize(String path, Path workDir) {
        if( !path.startsWith('/') )
            return path

        final dir = workDir.toString()
        if( !path.startsWith(dir) )
            throw new IllegalFileException("File `$path` is out of the scope of process working dir: $workDir")

        if( path.length()-dir.length()<2 )
            throw new IllegalFileException("Missing output file name")

        return path.substring(dir.size()+1)
    }

    @PackageScope
    String relativize(Path path, Path workDir) {
        if( !path.isAbsolute() )
            return glob ? FilePatternSplitter.GLOB.escape(path) : path

        if( !path.startsWith(workDir) )
            throw new IllegalFileException("File `$path` is out of the scope of process working dir: $workDir")

        if( path.nameCount == workDir.nameCount )
            throw new IllegalFileException("Missing output file name")

        final rel = path.subpath(workDir.getNameCount(), path.getNameCount())
        return glob ? FilePatternSplitter.GLOB.escape(rel) : rel
    }

    /**
     * Override the default to allow null as a value name
     * @return
     */
    String getName() {
        return nameObj ? super.getName() : null
    }

    @Override
    FileOutParam setPathQualifier(boolean flag) {
        pathQualifier = flag
        separatorChar = null
        return this
    }

    @Override
    boolean isPathQualifier() { pathQualifier }

    @Override
    FileOutParam setOptions(Map<String,?> opts) {
        (FileOutParam)super.setOptions(opts)
    }

}
