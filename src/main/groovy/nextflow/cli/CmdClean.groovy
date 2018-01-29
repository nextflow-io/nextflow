/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
import picocli.CommandLine

/**
 * Implements cache clean up command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
//@Parameters(commandDescription = "Clean up project cache and work directories")
@CommandLine.Command (name = "Clean", description ="Clean up project cache and work directories")
class CmdClean extends CmdBase implements CacheBase {

    static final public NAME = 'clean'

    //@Parameter(names=['-q', '-quiet'], description = 'Do not print names of files removed', arity = 0)
    @CommandLine.Option(names=['-q', '--quiet'], description = 'Do not print names of files removed', arity = '0')
    boolean quiet

    //@Parameter(names=['-f', '-force'], description = 'Force clean command', arity = 0)
    @CommandLine.Option(names=['-f', '--force'], description = 'Force clean command', arity = '0')
    boolean force

    //@Parameter(names=['-n', '-dry-run'], description = 'Print names of file to be removed without deleting them' , arity = 0)
    @CommandLine.Option(names=['-n', '--dry-run'], description = 'Print names of file to be removed without deleting them' , arity = '0')
    boolean dryRun

    //@Parameter(names='-after', description = 'Clean up runs executed after the specified one')
    @CommandLine.Option(names='--after', description = 'Clean up runs executed after the specified one')
    String after

    //@Parameter(names='-before', description = 'Clean up runs executed before the specified one')
    @CommandLine.Option(names='--before', description = 'Clean up runs executed before the specified one')
    String before

    //@Parameter(names='-but', description = 'Clean up all runs except the specified one')
    @CommandLine.Option(names='--but', description = 'Clean up all runs except the specified one')
    String but

    //@Parameter
    @CommandLine.Parameters(description = "")    //TODO is it mandatory?
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
        if( dryRun ) return

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
                println "Would remove ${record.workDir}"
            return
        }

        // decrement the ref count in the db
        def deleted = currentCacheDb.removeTaskEntry(hash)
        if( deleted ) {
            // delete folder
            if( deleteFolder(FileHelper.asPath(record.workDir))) {
                if(!quiet) println "Removed ${record.workDir}"
            }
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
    private boolean deleteFolder( Path folder ) {

        def result = true
        Files.walkFileTree(folder, new FileVisitor<Path>() {

            @Override
            FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) throws IOException {
                result ? FileVisitResult.CONTINUE : FileVisitResult.TERMINATE
            }

            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                if( !file.delete() ) {
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
                if( result && !dir.delete() ) {
                    result = false
                    if(!quiet) System.err.println "Failed to remove $dir"
                }

                result ? FileVisitResult.CONTINUE : FileVisitResult.TERMINATE
            }
        })

        return result
    }

}
