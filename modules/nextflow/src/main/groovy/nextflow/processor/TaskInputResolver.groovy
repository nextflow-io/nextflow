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

package nextflow.processor

import java.nio.file.Path
import java.util.regex.Matcher
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Nextflow
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.Executor
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.file.FilePorter
import nextflow.file.LogicalDataPath
import nextflow.script.ScriptType
import nextflow.script.params.FileInParam
import nextflow.script.params.v2.ProcessFileInput
import nextflow.util.ArrayBag
import nextflow.util.BlankSeparatedList
/**
 * Implements the resolution of input files for a task.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskInputResolver {

    private TaskRun task

    private FilePorter.Batch foreignFiles

    private Executor executor

    private int count = 0

    TaskInputResolver(TaskRun task, FilePorter.Batch foreignFiles, Executor executor) {
        this.task = task
        this.foreignFiles = foreignFiles
        this.executor = executor
    }

    /**
     * Resolve a file input for a legacy process.
     *
     * @param param
     * @param value
     */
    List<FileHolder> resolve(FileInParam param, Object value) {
        final ctx = task.context
        final normalized = normalizeInputToFiles(value, count, param.isPathQualifier())
        final resolved = expandWildcards( param.getFilePattern(ctx), normalized )

        // add to context if the path was declared with a variable name
        if( param.name )
            ctx.put( param.name, singleItemOrList(resolved, param.isSingle(), task.type) )

        count += resolved.size()

        return resolved
    }

    /**
     * Resolve a file input for a typed process.
     *
     * @param fileInput
     * @param value
     */
    List<FileHolder> resolve(ProcessFileInput fileInput, Object value) {
        final ctx = task.context
        final normalized = normalizeInputToFiles(value, count, true)
        final resolved = expandWildcards( fileInput.getFilePattern(ctx), normalized )

        count += resolved.size()

        return resolved
    }

    /**
     * Transform a value (i.e. path, collection, or map) by
     * replacing any source paths with staged paths.
     *
     * @param value
     * @param holders
     */
    Object normalizeValue(Object value, Map<Path,FileHolder> holders) {
        if( value instanceof Path ) {
            return normalizePath(value, holders)
        }

        if( value instanceof Collection ) {
            return value.collect { el -> normalizeValue(el, holders) }
        }

        if( value instanceof Map ) {
            return value.collectEntries { k, v -> [k, normalizeValue(v, holders)] }
        }

        return value
    }

    private Path normalizePath(Path value, Map<Path,FileHolder> holders) {
        return holders.containsKey(value)
            ? new TaskPath(holders[value])
            : value
    }

    protected List<FileHolder> normalizeInputToFiles( Object obj, int count, boolean coerceToPath ) {

        Collection allItems = obj instanceof Collection ? obj : [obj]
        def len = allItems.size()

        // use a bag so that cache hash key is not affected by file entries order
        def files = new ArrayBag<FileHolder>(len)
        for( def item : allItems ) {

            if( item instanceof Path || coerceToPath ) {
                final path = normalizeToPath(item)
                final target = executor.isForeignFile(path) ? foreignFiles.addToForeign(path) : path
                final holder = new FileHolder(path, target)
                files << holder
            }
            else {
                files << normalizeInputToFile(item, "input.${++count}")
            }
        }

        return files
    }

    protected static Path normalizeToPath( obj ) {
        if( obj instanceof LogicalDataPath )
            return obj.toTargetPath()

        if( obj instanceof Path )
            return obj

        if( obj == null )
            throw new ProcessUnrecoverableException("Path value cannot be null")
        
        if( !(obj instanceof CharSequence) )
            throw new ProcessUnrecoverableException("Not a valid path value type: ${obj.getClass().getName()} ($obj)")

        def str = obj.toString().trim()
        if( str.contains('\n') )
            throw new ProcessUnrecoverableException("Path value cannot contain a new-line character: $str")
        if( str.startsWith('/') )
            return FileHelper.asPath(str)
        if( FileHelper.getUrlProtocol(str) )
            return FileHelper.asPath(str)
        if( !str )
            throw new ProcessUnrecoverableException("Path value cannot be empty")
        
        throw new ProcessUnrecoverableException("Not a valid path value: '$str'")
    }

    /**
     * An input file parameter can be provided with any value other than a file.
     * This function normalize a generic value to a {@code Path} create a temporary file
     * in the for it.
     *
     * @param input The input value
     * @param altName The name to be used when a temporary file is created.
     * @return The {@code Path} that will be staged in the task working folder
     */
    protected static FileHolder normalizeInputToFile( Object input, String altName ) {
        /*
         * when it is a local file, just return a reference holder to it
         */
        if( input instanceof Path ) {
            return new FileHolder(input)
        }

        /*
         * default case, convert the input object to a string and save
         * to a local file
         */
        def source = input?.toString() ?: ''
        def result = Nextflow.tempFile(altName)
        result.text = source
        return new FileHolder(source, result)
    }

    /**
     * An input file name may contain wildcards characters which have to be handled coherently
     * given the number of files specified.
     *
     * @param name A file name with may contain a wildcard character star {@code *} or question mark {@code ?}.
     *  Only one occurrence can be specified for star or question mark wildcards.
     *
     * @param value Any value that have to be managed as an input files. Values other than {@code Path} are converted
     * to a string value, using the {@code #toString} method and saved in the local file-system. Value of type {@code Collection}
     * are expanded to multiple values accordingly.
     *
     * @return
     */
    protected static List<FileHolder> expandWildcards( String name, List<FileHolder> files ) {
        assert files != null

        // use an unordered so that cache hash key is not affected by file entries order
        final result = new ArrayBag(files.size())
        if( files.size()==0 ) { return result }

        if( !name || name == '*' ) {
            result.addAll(files)
            return result
        }

        if( !name.contains('*') && !name.contains('?') && files.size()>1 ) {
            /*
             * When name do not contain any wildcards *BUT* multiple files are provide
             * it is managed like having a 'star' at the end of the file name
             */
            name += '*'
        }

        for( int i=0; i<files.size(); i++ ) {
            def holder = files[i]
            def newName = expandWildcards0(name, holder.stageName, i+1, files.size())
            result << holder.withName( newName )
        }

        return result
    }

    protected static String expandWildcards0( String path, String stageName, int index, int size ) {

        String name
        String parent
        int p = path.lastIndexOf('/')
        if( p == -1 ) {
            parent = null
            name = path
        }
        else {
            parent = path.substring(0,p)
            name = path.substring(p+1)
        }

        if( name == '*' || !name ) {
            name = stageName
        }
        else {
            final stripWildcard = size<=1 // <-- string the start wildcard instead of expanding to an index number when the collection contain only one file
            name = replaceStarWildcards(name, index, stripWildcard)
            name = replaceQuestionMarkWildcards(name, index)
        }

        if( parent ) {
            parent = replaceStarWildcards(parent, index)
            parent = replaceQuestionMarkWildcards(parent, index)
            return "$parent/$name"
        }
        else {
            return name
        }
    }

    protected static String replaceStarWildcards(String name, int index, boolean strip=false) {
        name.replaceAll(/\*/, strip ? '' : String.valueOf(index))
    }

    private static final Pattern QUESTION_MARK = ~/(\?+)/

    protected static String replaceQuestionMarkWildcards(String name, int index) {
        def result = new StringBuffer()

        Matcher m = QUESTION_MARK.matcher(name)
        while( m.find() ) {
            def match = m.group(1)
            def repString = String.valueOf(index).padLeft(match.size(), '0')
            m.appendReplacement(result, repString)
        }
        m.appendTail(result)
        result.toString()
    }

    protected static Object singleItemOrList( List<FileHolder> items, boolean single, ScriptType type ) {
        assert items != null

        if( items.size() == 1 && single ) {
            return makePath(items[0],type)
        }

        def result = new ArrayList(items.size())
        for( int i=0; i<items.size(); i++ ) {
            result.add( makePath(items[i],type) )
        }
        return new BlankSeparatedList(result)
    }

    private static Path makePath( FileHolder holder, ScriptType type ) {
        if( type == ScriptType.SCRIPTLET ) {
            return new TaskPath(holder)
        }
        if( type == ScriptType.GROOVY) {
            // the real path for the native task needs to be fixed -- see #378
            return Path.of(holder.stageName)
        }
        throw new IllegalStateException("Unknown task type: $type")
    }
}
