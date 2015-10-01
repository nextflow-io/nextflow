/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
package nextflow.executor
import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.exception.ProcessException
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.processor.TaskRun
import nextflow.util.KryoHelper
import nextflow.util.RemoteSession
import org.apache.ignite.Ignite
import org.apache.ignite.IgniteCache
import org.apache.ignite.IgniteLogger
import org.apache.ignite.compute.ComputeJob
import org.apache.ignite.lang.IgniteCallable
import org.apache.ignite.resources.IgniteInstanceResource
import org.apache.ignite.resources.LoggerResource

/**
 * Models a task executed remotely in a Ignite cluster node
 *
 * @param < T > The type of the value returned by the {@code #call} method
 */
@CompileStatic
abstract class IgBaseTask<T> implements IgniteCallable<T>, ComputeJob {

    static final Map<UUID,GroovyClassLoader> classLoaderCache = new HashMap()

    @LoggerResource
    private transient IgniteLogger log

    @IgniteInstanceResource
    private transient Ignite grid

    /**
     * The client session identifier, it is required in order to access to
     * remote class-path
     */
    UUID sessionId

    /**
     * This field is used to transport the class attributes as a unique serialized byte array
     */
    private byte[] payload

    /**
     * Holds the class attributes in this map. Note: is defined as 'transient' because
     * the map content is serialized as a byte[] and saved to the {@code payload} field
     */
    private transient Map<String,Object> attrs = [:]

    /**
     * The local scratch dir where the task is actually executed in the remote node.
     * Note: is declared transient because it is valid only on the remote-side,
     * thus it do not need to be transported
     *
     */
    protected transient Path scratchDir

    /**
     * A temporary where all files are cached. The folder is deleted during instance shut-down
     */
    private static final Path localCacheDir = FileHelper.createLocalDir()

    protected static getLocalCacheDir() { localCacheDir }

    static {
        Runtime.getRuntime().addShutdownHook { localCacheDir.deleteDir() }
    }

    /**
     * Initialize the grid gain task wrapper
     *
     * @param task The task instance to be executed
     * @param The session unique identifier
     */
    protected IgBaseTask( TaskRun task, UUID sessionId ) {

        this.sessionId = sessionId

        attrs.taskId = task.id
        attrs.name = task.name
        attrs.workDir = task.workDir
        attrs.targetDir = task.targetDir
        attrs.inputFiles = [:]
        attrs.outputFiles = []

        // -- The a mapping of input files and target names
        def allFiles = task.getInputFiles().values()
        for( List<FileHolder> entry : allFiles ) {
            if( entry ) for( FileHolder it : entry ) {
                attrs.inputFiles[ it.stageName ] = it.storePath
            }
        }

        // -- the list of expected file names in the scratch dir
        attrs.outputFiles = task.getOutputFilesNames()

        payload = KryoHelper.serialize(attrs)
    }

    /** ONLY FOR TESTING PURPOSE */
    protected IgBaseTask() {}

    /**
     * @return The task unique ID
     */
    protected Object getTaskId() { attrs.taskId }

    /**
     * @return The task descriptive name (only for debugging)
     */
    protected String getName() { attrs.name }

    /**
     * @return The path where result files have to be copied
     */
    protected Path getTargetDir() { (Path)attrs.targetDir }

    /**
     * @return The task working directory i.e. the folder containing the scripts files, but
     * it is not the actual task execution directory
     */
    protected Path getWorkDir() { (Path)attrs.workDir }

    /**
     * @return The a mapping of input files and target names
     */
    Map<String,Path> getInputFiles() { (Map<String,Path>)attrs.inputFiles }

    /**
     * @return the list of expected file names in the scratch dir
     */
    List<String> getOutputFiles() { (List<String>)attrs.outputFiles }

    /**
     * Copies to the task input files to the execution folder, that is {@code scratchDir}
     * folder created when this method is invoked
     *
     */
    protected void stage() {

        if( attrs == null && payload )
            attrs = (Map<String,Object>)KryoHelper.deserialize(payload)

        // create a local scratch dir
        scratchDir = FileHelper.createLocalDir()

        if( !inputFiles )
            return

        // move the input files there
        for( Map.Entry<String,Path> entry : inputFiles.entrySet() ) {
            final fileName = entry.key
            final source = entry.value
            final cached = FileHelper.getLocalCachePath(source,localCacheDir, sessionId)
            final staged = scratchDir.resolve(fileName)
            log?.debug "Task ${getName()} > staging path: '${source}' to: '$staged'"
            Files.createSymbolicLink(staged, cached)
        }
    }


    /**
     * Copy back the task output files from the execution directory in the local node storage
     * to the task {@code targetDir}
     */
    protected void unstage() {
        log?.debug "Unstaging file names: $outputFiles"

        if( !outputFiles )
            return

        // create a bash script that will copy the out file to the working directory
        if( !Files.exists(targetDir) )
            Files.createDirectories(targetDir)

        for( String name : outputFiles ) {
            try {
                copyToTargetDir(name, scratchDir, targetDir)
            }
            catch( IOException e ) {
                log.error("Unable to copy result file: $name to target dir", e)
            }
        }
    }

    /**
     * Copy the file with the specified name from the task execution folder
     * to the {@code targetDir}
     *
     * @param filePattern A file name relative to the {@code scratchDir}.
     *        It can contain globs wildcards
     */
    protected void copyToTargetDir( String filePattern, Path from, Path to ) {

        def type = filePattern.contains('**') ? 'file' : 'any'

        FileHelper.visitFiles( from, filePattern, type: type ) { Path it ->
            final rel = from.relativize(it)
            it.copyTo(to.resolve(rel))
        }
    }


    /**
     * Invoke the task execution. It calls the following methods in this sequence: {@code stage}, {@code execute0} and {@code unstage}
     *
     * @return The {@code execute0} result value
     * @throws nextflow.exception.ProcessException
     */
    @Override
    final T call() throws Exception {
        try {
            /*
             * stage the input files in the working are`
             */
            stage()

            /*
             * execute the task
             */
            final T result = execute0()

            /*
             * copy back the result files to the shared area
             */
            unstage()

            // return the exit status eventually
            return result
        }
        catch( Exception e ) {
            log.error("Cannot execute task > $name", e)
            throw new ProcessException(e)
        }

    }

    /**
     * Just a synonym for {@code #call}
     *
     * @return The value returned by the task execution
     */
    final Object execute() {
        call()
    }

    /**
     * The actual task executor code provided by the extending subclass
     *
     * @return The value returned by the task execution
     */
    protected abstract T execute0()

    /**
     * Lookup the {@link RemoteSession} object for the given session ID
     *
     * @param sessionId The remote session ID
     * @param grid The current {@link Ignite} instance
     * @return The associated {@link RemoteSession} object for the specified session ID
     * @throws IllegalStateException when no session is found for the specified session ID
     */
    protected RemoteSession getSessionFor( UUID sessionId ) {
        assert sessionId
        IgniteCache<UUID, RemoteSession> allSessions = grid.cache( IgGridFactory.SESSIONS_CACHE )

        if( !allSessions )
            throw new IllegalStateException('Missing session cache object')

        def session = allSessions.get(sessionId)
        if( !session )
            throw new IllegalStateException("Missing session object for id: $sessionId")

        return session
    }

    /**
     * Create a {@link ClassLoader} object for the specified session ID
     *
     * @param sessionId
     * @param grid
     * @return
     */
    protected ClassLoader getClassLoaderFor( UUID sessionId ) {
        assert sessionId

        classLoaderCache.getOrCreate(sessionId) {

            final allSessions = (IgniteCache<UUID, RemoteSession>)grid.cache( IgGridFactory.SESSIONS_CACHE )
            if( !allSessions )
                throw new IllegalStateException('Missing session cache object')

            final session = allSessions.get(sessionId)
            if( !session )
                throw new IllegalStateException("Missing session object for id: $sessionId")

            final result = new GroovyClassLoader()
            session.classpath.each { Path file ->
                log.debug "Adding to classpath: $file"
                result.addClasspath(file.toAbsolutePath().toString())
            }

            return result
        }
    }
}