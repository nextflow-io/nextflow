/*
 * Copyright 2013-2026, Seqera Labs
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
 */

package io.seqera.wave.plugin

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.trace.TraceObserverV2
import nextflow.trace.TraceObserverFactoryV2
/**
 * Factory class for wave session
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WaveFactory implements TraceObserverFactoryV2 {

    static final String SEQERA_EXECUTOR = 'seqera'

    @Override
    Collection<TraceObserverV2> create(Session session) {
        shouldEnable(session)
        return Collections.<TraceObserverV2>emptyList()
    }

    @Memoized // <-- declare as memoized to make it's invoked only once
    static boolean shouldEnable(Session session) {
        final config = session.config
        final wave = (Map)config.wave ?: new HashMap<>(1)
        final fusion = (Map)config.fusion ?: new HashMap<>(1)

        if( SysEnv.get('NXF_DISABLE_WAVE_SERVICE') ) {
            log.debug "Detected NXF_DISABLE_WAVE_SERVICE environment variable - Turning off Wave service"
            wave.enabled = false
            return false
        }

        // the Seqera executor hosts provide the Fusion client, therefore Wave is optional
        if( fusion.enabled && isSeqeraExecutor(config) ) {
            if( wave.enabled )
                enableBundleProjectResources(session, wave, 'Fusion')
        }
        else if( fusion.enabled ) {
            checkWaveRequirement(session, wave, 'Fusion')
        }
        if( isAwsBatchFargateMode(config) ) {
            checkWaveRequirement(session, wave, 'Fargate')
        }
        return wave.enabled==true
    }

    static private void checkWaveRequirement(Session session, Map wave, String feature) {
        if( !wave.enabled ) {
            throw new AbortOperationException("$feature feature requires enabling Wave service")
        }
        else {
            enableBundleProjectResources(session, wave, feature)
        }
    }

    static private void enableBundleProjectResources(Session session, Map wave, String feature) {
        log.debug "Detected $feature enabled -- Enabling bundle project resources -- Disabling upload of remote bin directory"
        wave.bundleProjectResources = true
        session.disableRemoteBinDir = true
    }

    /**
     * Whether the run uses the Seqera executor, i.e. it is the default executor for every process.
     * A per-process ({@code withName}/{@code withLabel}) executor is not considered.
     *
     * @param config the session config
     * @return {@code true} when the default executor is {@code seqera}
     */
    static boolean isSeqeraExecutor(Map config) {
        final executor = config.navigate('process.executor') ?: config.navigate('executor.name')
        return SEQERA_EXECUTOR == executor?.toString()
    }

    static boolean isAwsBatchFargateMode(Map config) {
        return 'fargate'.equalsIgnoreCase(config.navigate('aws.batch.platformType') as String)
    }
}
