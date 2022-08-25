package io.seqera.wave.plugin

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.AbortRunException
import nextflow.trace.TraceObserver
import nextflow.trace.TraceObserverFactory
/**
 * Factory class for wave session
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WaveFactory implements TraceObserverFactory {

    @Override
    Collection<TraceObserver> create(Session session) {
        final config = session.config
        final wave = (Map)config.wave
        final fusion = (Map)config.fusion
        if( fusion?.enabled ) {
            if( !wave?.enabled ) throw new AbortRunException("Fusion feature requires enabling Wave service")
            log.debug "Detected Fusion enabled -- Enabling bundle project resources -- Disabling upload of remote bin directory"
            wave.bundleProjectResources = true
            config.disableRemoteBinDir = true
        }
        return Collections.emptyList()
    }
}
