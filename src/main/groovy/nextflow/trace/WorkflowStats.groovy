package nextflow.trace

import java.text.DecimalFormat

/**
 * Value object representing the workflow execution statistics
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WorkflowStats {

    /**
     * Overall workflow compute time (CPUs-seconds) for task executed successfully
     */
    long timeSucceed

    /**
     * Overall compute time (CPUs-seconds) for cached tasks
     */
    long timeCached

    /**
     * Overall compute time (CPUs-seconds) for failed tasks
     */
    long timeFailed

    /**
     * Task successfully completed
     */
    long completed

    /**
     * Task cached
     */
    long cached

    /**
     * Task terminated with an error
     */
    long failed

    /**
     * Task executions terminated with an error which was ignored
     */
    long ignored

    /**
     * @return A formatted string representing the overall execution time as CPU-Hours
     */
    protected String getComputeTime() {

        def fmt = new DecimalFormat("0.#")

        def total = (timeSucceed + timeCached + timeFailed)
        def result = String.format('%.1f', total/3600)
        if( timeCached || timeFailed ) {
            result += ' ('
            def items = []
            if( timeCached )
                items << fmt.format(timeCached/total*100) + '% cached'
            if( timeFailed )
                items << fmt.format(timeFailed/total*100) + '% failed'
            result += items.join(', ')
            result += ')'
        }

        return result
    }
}
