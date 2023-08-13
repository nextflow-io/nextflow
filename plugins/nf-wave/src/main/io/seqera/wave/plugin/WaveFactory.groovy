/*
 * Copyright 2013-2023, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

package io.seqera.wave.plugin

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
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
        final wave = (Map)config.wave ?: new HashMap<>(1)
        final fusion = (Map)config.fusion ?: new HashMap<>(1)

        if( SysEnv.get('NXF_DISABLE_WAVE_SERVICE') ) {
            log.debug "Detected NXF_DISABLE_WAVE_SERVICE environment variable - Turning off Wave service"
            wave.enabled = false
            return List.of()
        }
        
        if( fusion.enabled ) {
            if( !wave.enabled ) {
                throw new AbortOperationException("Fusion feature requires enabling Wave service")
            }
            else {
                log.debug "Detected Fusion enabled -- Enabling bundle project resources -- Disabling upload of remote bin directory"
                wave.bundleProjectResources = true
                session.disableRemoteBinDir = true
            }
        }

        final observer = new WaveObserver(session)
        return wave.enabled && observer.reportOpts().enabled()
                ? List.<TraceObserver>of(observer)
                : List.<TraceObserver>of()
    }
}
