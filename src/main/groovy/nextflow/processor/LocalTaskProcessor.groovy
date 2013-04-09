/*
 * Copyright (c) 2012, the authors.
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

package nextflow.processor
import groovy.io.FileType
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import nextflow.util.ByteDumper
import nextflow.util.Duration
import org.apache.commons.io.IOUtils
import org.codehaus.groovy.runtime.IOGroovyMethods

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
class LocalTaskProcessor extends AbstractTaskProcessor {

    private Duration maxDuration

    /**
     * The max duration time allowed for the job to be executed, this value sets the '-l h_rt' squb command line option.
     *
     * @param duration0 The max allowed time expressed as duration string, Accepted units are 'min', 'hour', 'day'.
     *                  For example {@code maxDuration '30 min'}, {@code maxDuration '10 hour'}, {@code maxDuration '2 day'}
     */
    LocalTaskProcessor maxDuration( String duration0 ) {
        this.maxDuration = new Duration(duration0)
        return this
    }


    /**
     * Run a system executable script
     *
     * @param script
     * @return
     */
    protected void launchTask( TaskDef task )  {
        assert task
        assert task.@script

        task.workDirectory = shareWorkDir ? new File(session.workDirectory, name) : new File(session.workDirectory, "$name-${task.index}")
        File scratch = task.workDirectory

        if ( !scratch.exists() && !scratch.mkdirs() ) {
            throw new IOException("Unable to create task work directory: '${scratch}'")
        }

        log.debug "Lauching task > ${task.name} -- scratch folder: $scratch"

        // -- create the command script file
        def scriptFile = new File(scratch, '.command.sh')
        scriptFile.createNewFile()
        scriptFile.text = task.script.toString().stripIndent()

        // -- save the reference to the scriptFile
        task.script = scriptFile

        ProcessBuilder builder = new ProcessBuilder()
                .directory(scratch)
                .command( shell, scriptFile.toString() )
                .redirectErrorStream(true)

        // -- configure the job environment
        builder.environment().putAll(getProcessEnvironment())

        // -- start the execution and notify the event to the monitor
        Process process = builder.start()
        task.status = TaskDef.Status.RUNNING

        // -- copy the input value to the process standard input
        if( task.input != null ) {
            pipeInput( task, process )
        }

        File fileOut = new File(scratch, '.command.out')
        ByteDumper dumper = null
        try {
            // -- print the process out if it is not capture by the output
            //    * The byte dumper uses a separate thread to capture the process stdout
            //    * The process stdout is captured in two condition:
            //      when the flag 'echo' is set or when it goes in the output channel (outputs['-'])
            //
            BufferedOutputStream streamOut = new BufferedOutputStream( new FileOutputStream(fileOut) )

            def handler = { byte[] data, int len ->
                streamOut.write(data,0,len)
                if( echo ) System.out.print(new String(data,0,len))
            }
            dumper = new ByteDumper(process.getInputStream(), handler)
            dumper.setName("dumper-$name")
            dumper.start()

            // -- wait the the process completes
            if( maxDuration ) {
                log.debug "Running task > ${task.name} -- waiting max: ${maxDuration}"
                process.waitForOrKill(maxDuration.toMillis())
                task.exitCode = process.exitValue()
            }
            else {
                log.debug "Running task > ${task.name} -- wait forever"
                task.exitCode = process.waitFor()
            }

            log.debug "Task completed > ${task.name} -- exit code: ${task.exitCode}; success: ${task.exitCode in validExitCodes}"

            dumper?.await(500)
            streamOut.close()

        }
        finally {
            dumper?.terminate()
            task.workDirectory = scratch
            task.output = fileOut

            IOUtils.closeQuietly(process.in)
            IOUtils.closeQuietly(process.out)
            IOUtils.closeQuietly(process.err)
            process.destroy()
        }

    }



    @Override
    protected List<File> collectResultFile( File scratchPath, String name ) {
        assert scratchPath
        assert name

        // replace any wildcards characters
        // TODO give a try to http://code.google.com/p/wildcard/  -or- http://commons.apache.org/io/
        String filePattern = name.replace("?", ".?").replace("*", ".*?")

        if( filePattern == name ) {
            // TODO check that the file exists (?)
            return [ new File(scratchPath,name) ]
        }

        // scan to find the file with that name
        List files = []
        scratchPath.eachFileMatch(FileType.FILES, ~/$filePattern/ ) { File it -> files << it}

        // TODO ++ what if expected files are missing?
        return files
    }

    /**
     * Pipe the {@code TaskDef#input} to the {@code Process}
     *
     * @param task The current task to be executed
     * @param process The system process that will run the task
     */
    protected void pipeInput( TaskDef task, Process process ) {

        Thread.start {
            try {
                IOGroovyMethods.withStream(new BufferedOutputStream(process.getOutputStream())) {writer -> writer << task.input}
            }
            catch( Exception e ) {
                log.warn "Unable to pipe input data for task: ${task.name}"
            }
        }


    }




}
