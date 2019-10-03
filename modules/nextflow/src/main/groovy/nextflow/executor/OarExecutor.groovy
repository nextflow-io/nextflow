// Copyright?

package nextflow.executor
import java.nio.file.Path
import java.util.regex.Pattern

import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
/**
 * Processor for OAR resource manager (DRAFT)
 *
 * See https://oar.imag.fr/
 *
 *
 * @author Maxime Vallée <maxime.vallee@chu-lyon.fr>
 */
@Slf4j
class OarExecutor extends AbstractGridExecutor {

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << '-d' << quote(task.workDir)
        result << '-n' << getJobNameFor(task)
        result << '-O' << quote(task.workDir.resolve(TaskRun.CMD_OUTFILE))     
        result << '-E' << quote(task.workDir.resolve(TaskRun.CMD_ERRFILE))     

		// Might be difficult, but all parameters are passed in one argument, like :
		// 1 core on 2 nodes on the same cluster with 16384 MB of memory and Infiniband 20G + 1 cpu on 2 nodes on the same switch with 8 cores processors for a walltime of 4 hours
		// oarsub -I -l "{memnode=16384 and ib_rate='20'}/cluster=1/nodes=2/core=1+{cpucore=8}/switch=1/nodes=2/cpu=1,walltime=4:0:0"
		// Warning
		// walltime must always be the last argument of -l <...>

        if( task.config.cpus > 1 ) {
            result << '-l cpu=' << task.config.cpus.toString() // NOT SURE AT ALL
        }

        if( task.config.getMemory() ) {
            result << '-l memnode=' << task.config.getMemory().toMega().toString() // NOT SURE AT ALL
        }

        if( task.config.time ) {
            result << '-l walltime=' << task.config.getTime().format('HH:mm:ss') // NOT SURE AT ALL
        }

        // the requested queue name
        if( task.config.queue ) {
            result << '-q' << (task.config.queue.toString())
        }

        // -- at the end append the command script wrapped file name
        if( task.config.clusterOptions ) {
            result << task.config.clusterOptions.toString() << ''
        }

        return result
    }

    String getHeaderToken() { '#OAR' }

    /**
     * The command line to submit this job
     *
     * @param task The {@link TaskRun} instance to submit for execution to the cluster
     * @param scriptFile The file containing the job launcher script
     * @return A list representing the submit command line
     */
    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {
		// Scripts need to be executable
		scriptFile.setPermissions(7,0,0)
        return ["oarsub", "-S", "./${scriptFile.getName()}"]
    }

    /**
     * Parse the string returned by the {@code oarsub} command and extract the job ID string
     *
     * @param text The string returned when submitting the job
     * @return The actual job ID string
     */
    @Override
    def parseJobId(String text) {
    def pattern = ~ /OAR_JOB_ID=(\d+)/
    for( String line : text.readLines() ) {
        def m = pattern.matcher(line)
        if( m.matches() ) {
            return m.group(1).toString()
            }
        }
        throw new IllegalStateException("Invalid OAR submit response:\n$text\n\n")
    }

    @Override
    protected List<String> getKillCommand() { ['oardel'] }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        // To have a parsable list of jobs in queue by user
		// see page 21 http://oar.imag.fr/docs/2.5/OAR-Documentation.pdf
		String cmd = 'oarstat -f'
        if( queue ) cmd += ' ' + queue
        return ['sh','-c', "$cmd | egrep '(Job_Id:|state =)' ".toString()]
    }

    /*
     *  Maps OAR job status to nextflow status
     *  see page 134 http://oar.imag.fr/docs/2.5/OAR-Documentation.pdf
     */
    static private Map DECODE_STATUS = [
            'toLaunch': QueueStatus.RUNNING,			// (the OAR scheduler has attributed some nodes to the job. So it will be launched.)
            'Launching': QueueStatus.RUNNING,			// (OAR has launched the job and will execute the user command on the first node.)
            'Running': QueueStatus.RUNNING,				// (the user command is executing on the first node.)
            'Finishing': QueueStatus.RUNNING,			// (the user command has terminated and OAR is doing work internally)
            'Waiting': QueueStatus.PENDING,				// (the job is waiting OAR scheduler decision.)
            'toAckReservation': QueueStatus.PENDING,	// (the OAR scheduler must say “YES” or “NO” to the waiting oarsub command because it requested a reservation.)
            'Hold': QueueStatus.HOLD,					// (user or administrator wants to hold the job (oarhold command). So it will not be scheduled by the system.)
            'Suspended': QueueStatus.HOLD,				// (the job was in Running state and there was a request (oarhold with “-r” option) to suspend this job. In this state other jobs can be scheduled on the same resources (these resources has the “suspended_jobs” field to “YES”).)
            'Error': QueueStatus.ERROR,					// (a problem has occurred.)
            'toError': QueueStatus.ERROR,				// (something wrong occurred and the job is going into the error state.)
            'Terminated': QueueStatus.DONE,				// (the job has terminated normally.)
    ]

    protected QueueStatus decode(String status) {
        DECODE_STATUS.get(status)
    }

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        final JOB_ID = 'Job_Id:'
        final JOB_STATUS = 'state ='
        final result = [:]

        String id = null
        String status = null
        text.eachLine { line ->
            if( line.startsWith(JOB_ID) ) {
                id = fetchValue(JOB_ID, line)
            }
            else if( id ) {
                status = fetchValue(JOB_STATUS, line)
            }
            result.put( id, decode(status) ?: QueueStatus.UNKNOWN )
        }

        return result
    }

    static String fetchValue( String prefix, String line ) {
        final p = line.indexOf(prefix)
        return p!=-1 ? line.substring(p+prefix.size()).trim() : null
    }


}
