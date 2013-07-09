package nextflow.executor
import groovy.transform.InheritConstructors
import nextflow.processor.TaskRun
/**
 * Processor for SLURM resource manager (DRAFT)
 *
 * See http://computing.llnl.gov/linux/slurm/
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class SlurmExecutor extends AbstractGridExecutor {


    /**
     * -c, --cpus-per-task=ncpus   number of cpus required per task
     * -D, --chdir=path            change remote current working directory
     * -e, --error=err             location of stderr redirection
     * -E, --preserve-env          env vars for node and task counts override
     * -i, --input=in              location of stdin redirection
     * -J, --job-name=jobname      name of job
     * -o, --output=out            location of stdout redirection
     * -Q, --quiet                 quiet mode (suppress informational messages)
     * -t, --time=minutes          time limit
     * -u, --unbuffered            do not line-buffer stdout/err
     * --mem=MB                minimum amount of real memory
     * --mincpus=n             minimum number of logical processors (threads) per node
     * --tmp=MB                minimum amount of temporary disk
     * --mem-per-cpu=MB        maximum amount of real memory per allocated cpu required by the job. --mem >= --mem-per-cpu if --mem is specified.
     *
     * @param task
     * @return
     */
    @Override
    protected List<String> getSubmitCommandLine(TaskRun task) {

        final result = new ArrayList<String>()

        result << 'srun'
        result << '-D' << task.workDirectory
        result << '-J' << "nf-${task.processor.name}-${task.index}"
        result << '-E'

        if( taskConfig.maxDuration ) {
            result << '-t' << taskConfig.maxDuration.format('HH:mm:ss')
        }

        // -- at the end append the command script wrapped file name
        if( taskConfig.gridNativeOptions ) {
            result.addAll( getGridNativeOptionsAsList() )
        }

        // -- last entry to 'script' file name
        // replace with the 'shell' attribute
        result << 'bash' << JOB_SCRIPT_FILENAME

    }

//
//    @Override
//    protected String getJobTempFolder() {
//        '[ ! -z $TMPDIR ] && TMPDIR="$TMPDIR/`tr -dc A-Za-z0-9_ < /dev/urandom | head -c 9`" && mkdir -p $TMPDIR && cd $TMPDIR\npwd'
//    }


}
