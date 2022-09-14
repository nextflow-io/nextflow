/*
 * Copyright 2020-2022, Seqera Labs
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
        final wave = (Map)config.wave
        final fusion = (Map)config.fusion
        if( fusion?.enabled ) {
            if( !wave?.enabled ) throw new AbortOperationException("Fusion feature requires enabling Wave service")
            log.debug "Detected Fusion enabled -- Enabling bundle project resources -- Disabling upload of remote bin directory"
            wave.bundleProjectResources = true
            session.disableRemoteBinDir = true
        }
        return Collections.emptyList()
    }
}
