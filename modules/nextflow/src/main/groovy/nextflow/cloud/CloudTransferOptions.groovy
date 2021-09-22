package nextflow.cloud

import nextflow.util.Duration

/**
 * Define common cloud data transfer options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface CloudTransferOptions {

    static final public int MAX_TRANSFER = 4

    static final public int MAX_TRANSFER_ATTEMPTS = 1

    static final public Duration DEFAULT_DELAY_BETWEEN_ATTEMPTS = Duration.of('10sec')

    int getMaxParallelTransfers()
    int getMaxTransferAttempts()
    Duration getDelayBetweenAttempts()

}
