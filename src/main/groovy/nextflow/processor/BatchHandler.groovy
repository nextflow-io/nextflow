package nextflow.processor

/**
 * Defines the contract for {@link TaskHandler} classes that need
 * to aggregate multiple operation to optimise the interaction with
 * a remote execution system
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
trait BatchHandler<K,V> {

    /**
     * Defines a {@link BatchContext} for the this handler
     *
     * @param context The {@link BatchContext} instance
     */
    abstract void batch( BatchContext<K,V> context )

}
