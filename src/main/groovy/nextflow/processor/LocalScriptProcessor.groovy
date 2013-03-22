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
import java.util.concurrent.CountDownLatch
import java.util.concurrent.TimeUnit

import groovy.io.FileType
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
class LocalScriptProcessor extends AbstractScriptProcessor {

    /**
     * Run a system executable script
     *
     * @param script
     * @return
     */
    protected void runScript( def script, TaskDef task )  {
        assert script

        task.workDirectory = shareWorkDirectory ? new File(session.workDirectory, name) : new File(session.workDirectory, "$name-${task.index}")
        File scratch = task.workDirectory

        if ( !scratch.exists() && !scratch.mkdirs() ) {
            throw new IOException("Unable to create task work directory: '${scratch}'")
        }

        log.debug "Processor '$name' running script using scratch folder: $scratch"

        // -- create the command script file
        def scriptFile = new File(scratch, '.command.sh')
        scriptFile.createNewFile()
        scriptFile.text = script.toString().stripIndent()

        ProcessBuilder builder = new ProcessBuilder()
                .directory(scratch)
                .command( shell, scriptFile.toString() )
                .redirectErrorStream(true)

        // -- configure the job environment
        builder.environment().putAll(getProcessEnv())

        // -- start the execution and notify the event to the monitor
        Process process = builder.start()
        task.status = TaskDef.Status.RUNNING

        // -- print the process out if it is not capture by the output
        //    * The byte dumper uses a separate thread to capture the process stdout
        //    * The process stdout is captured in two condition:
        //      when the flag 'echo' is set or when it goes in the output channel (outputs['-'])
        //
        ByteDumper dumper = null
        def buffer = outputs.containsKey('-') ? new ByteArrayOutputStream() : null
        if ( echo || buffer ) {
            def handler = { byte[] data, int len ->
                if( echo ) System.out.print(new String(data,0,len))
                if( buffer ) buffer.write(data,0,len)
            }
            dumper = new ByteDumper(process.getInputStream(), handler)
            dumper.setName("dumper-$name")
            dumper.start()
        }


        try {

            // -- wait the the process completes
            //    and collect the result
            task.exitCode = process.waitFor()
            dumper?.await(500)


        }
        finally {
            dumper?.terminate()
            task.workDirectory = scratch
            task.output = buffer ? new String(buffer.toByteArray()): null

        }

    }

    protected Map<String,String> getProcessEnv() {
        def result = new HashMap( session.env )
        if( environment ) {
            result.putAll(environment)
        }

        return result
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


    private static class ByteDumper extends Thread {

        InputStream fInput
        boolean fTerminated
        CountDownLatch barrier = new CountDownLatch(1)
        Closure fCallback

        public ByteDumper(InputStream _in, Closure callback ) {
            assert _in
            assert callback

            this.fInput = new BufferedInputStream(_in);
            this.fCallback = callback
        }

        /**
         *  Interrupt the dumper thread
         */
        def void terminate() { fTerminated = true }

        /**
         * Await that the thread finished to read the process stdout
         *
         * @param millis Maximum time (in millis) to await
         */
        def void await(long millis=0) {
            if( millis ) {
                barrier.await(millis, TimeUnit.MILLISECONDS)
            }
            else {
                barrier.await()

            }
        }


        @Override
        public void run() {

            try {
                consume()
            }
            finally{
                barrier.countDown()
            }

        }

        public void consume() {
            byte[] buf = new byte[8192];
            int next;
            try {
                while ((next = fInput.read(buf)) != -1 && !fTerminated ) {
                    fCallback.call(buf, next)
                }
            } catch (IOException e) {
                throw new RuntimeException("exception while dumping process stream", e);
            }
        }
    }

}
