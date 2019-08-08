package nextflow.trace

import nextflow.Session

/**
 * Factory class creating {@link TraceObserver} instances
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface TraceObserverFactory {

    /**
     * Register the observer on the current session object
     *
     * @param session The current {@link nextflow.Session} instance
     * @return One or more instances of {@link TraceObserver} objects
     */
    Collection<TraceObserver> create(Session session)

}
