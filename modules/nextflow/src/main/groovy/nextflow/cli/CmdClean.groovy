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

package nextflow.cli
import java.nio.file.FileVisitResult
import java.nio.file.FileVisitor
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.CacheDB
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.trace.TraceRecord
import nextflow.util.HistoryFile.Record

/**
 * Implements cache clean up command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Lorenz Gerber <lorenzottogerber@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Clean up project cache and work directories")
class CmdClean extends CmdBase implements CacheBase {

    static final public NAME = 'clean'

    @Parameter(names=['-q', '-quiet'], description = 'Do not print names of files removed', arity = 0)
    boolean quiet

    @Parameter(names=['-f', '-force'], description = 'Force clean command', arity = 0)
    boolean force

    @Parameter(names=['-n', '-dry-run'], description = 'Print names of file to be removed without deleting them' , arity = 0)
    boolean dryRun

    @Parameter(names='-after', description = 'Clean up runs executed after the specified one')
    String after

    @Parameter(names='-before', description = 'Clean up runs executed before the specified one')
    String before

    @Parameter(names='-but', description = 'Clean up all runs except the specified one')
    String but

    @Parameter(names=['-k', '-keep-logs'], description = 'Removes only temporary files but retains execution log entries and metadata')
    boolean keepLogs

    @Parameter
    List<String> args

    private CacheDB currentCacheDb

    private Map<HashCode, Short> dryHash = new HashMap<>()

    /**
     * @return The name of this command {@code clean}
     */
    @Override
    String getName() {
        return NAME
    }

    /**
     * Command entry method
     */
    @Override
    void run() {
        init()
        validateOptions()

        listIds().each { entry -> cleanup(entry)}

    }

    /**
     * Extra CLI option validation
     */
    private void validateOptions() {

        if( !dryRun && !force )
            throw new AbortOperationException("Neither -f or -n specified -- refused to clean")
    }

    /**
     * Given a history entry clean up execution cache, deleting
     * task work directories and cache DB records
     *
     * @param entry
     *      A {@link Record} object representing a row in the history log file
     */
    private void cleanup(Record entry) {
        currentCacheDb = cacheFor(entry).openForRead()
        // -- remove each entry and work dir
        currentCacheDb.eachRecord(this.&removeRecord)
        // -- close the cache
        currentCacheDb.close()

        // -- STOP HERE !
        if( dryRun || keepLogs ) return

        // -- remove the index file
        currentCacheDb.deleteIndex()
        // -- remove the session from the history file
        history.deleteEntry(entry)
        // -- check if exists another history entry for the same session
        if( !history.checkExistsById(entry.sessionId)) {
            currentCacheDb.drop()
        }
    }

    /**
     * Check if a tasks can be removed during a dry-run simulation.
     *
     * @param hash
     *      The task unique hash code
     * @param refCount
     *      The number of times the task cache is references by other run instances
     * @return
     *      {@code true} when task will be removed by the clean command, {@code false} otherwise i.e.
     *      entry cannot be deleted because is referenced by other run instances
     */
    private boolean wouldRemove(HashCode hash, Integer refCount) {

        if( dryHash.containsKey(hash) ) {
            refCount = dryHash.get(hash)-1
        }

        if( refCount == 1 ) {
            dryHash.remove(hash)
            return true
        }
        else {
            dryHash.put(hash, (short)refCount)
            return false
        }

    }

    /**
     * Delete task cache entry
     *
     * @param hash The task unique hash code
     * @param record The task {@link TraceRecord}
     * @param refCount The number of times the task cache is references by other run instances
     */
    private void removeRecord(HashCode hash, TraceRecord record, int refCount) {
        if( dryRun ) {
            if( wouldRemove(hash,refCount) )
                printMessage(record.workDir,true)
            return
        }

        // decrement the ref count in the db
        def proceed = keepLogs || currentCacheDb.removeTaskEntry(hash)
        if( proceed ) {
            // delete folder
            if( deleteFolder(FileHelper.asPath(record.workDir), keepLogs)) {
                if(!quiet) printMessage(record.workDir,false)
            }

        }
    }

    private printMessage(String path, boolean dryRun) {
        if( dryRun ) {
            println keepLogs ? "Would remove temp files from ${path}" : "Would remove ${path}"
        }
        else {
            println keepLogs ? "Removed temp files from ${path}" : "Removed ${path}"
        }
    }

    /**
     * Traverse a directory structure and delete all the content
     *
     * @param folder
     *      The directory to delete
     * @return
     *      {@code true} in the directory was removed, {@code false}  otherwise
     */
    private boolean deleteFolder( Path folder, boolean keepLogs ) {

        def result = true
        Files.walkFileTree(folder, new FileVisitor<Path>() {

            @Override
            FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) throws IOException {
                result ? FileVisitResult.CONTINUE : FileVisitResult.TERMINATE
            }

            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {

                final canDelete = !keepLogs || ( keepLogs &&  !(file.name.startsWith('.command.')  || file.name == '.exitcode'))
                if( canDelete && !file.delete() ) {
                    result = false
                    if(!quiet) System.err.println "Failed to remove $file"
                }

                result ? FileVisitResult.CONTINUE : FileVisitResult.TERMINATE
            }

            @Override
            FileVisitResult visitFileFailed(Path file, IOException exc) throws IOException {
                FileVisitResult.CONTINUE
            }

            @Override
            FileVisitResult postVisitDirectory(Path dir, IOException exc) throws IOException {
                if( !keepLogs && result && !dir.delete() ) {
                    result = false
                    if(!quiet) System.err.println "Failed to remove $dir"
                }
                
                result ? FileVisitResult.CONTINUE : FileVisitResult.TERMINATE
            }
        })

        return result
    }

}
