package nextflow.util
import groovy.transform.CompileStatic
import groovyx.gpars.scheduler.Pool
import groovyx.gpars.util.PoolFactory
/**
 * A custom pool factory that will create instances of
 * {@link UnboundThreadPool}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CustomPoolFactory implements PoolFactory {

    @Override
    Pool createPool() {
        new UnboundThreadPool(true, retrievePoolSize())
    }

    @Override
    Pool createPool(boolean daemon) {
        new UnboundThreadPool(daemon, retrievePoolSize())
    }

    @Override
    Pool createPool(int numberOfThreads) {
        new UnboundThreadPool(true, numberOfThreads)
    }

    @Override
    Pool createPool(boolean daemon, int numberOfThreads) {
        new UnboundThreadPool(daemon, numberOfThreads)
    }

    static private int retrievePoolSize() {
        int size = Runtime.getRuntime().availableProcessors() * 8
        final String poolSizeValue = System.getProperty("gpars.poolsize")
        if( poolSizeValue && poolSizeValue.isInteger() ) {
            size = poolSizeValue.toInteger()
        }

        return size
    }
}
