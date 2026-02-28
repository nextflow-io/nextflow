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

package nextflow.trace

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.trace.config.DagConfig
import nextflow.trace.config.ReportConfig
import nextflow.trace.config.TimelineConfig
import nextflow.trace.config.TraceConfig

/**
 * Creates Nextflow observes object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DefaultObserverFactory implements TraceObserverFactoryV2 {

    private Session session

    @Override
    Collection<TraceObserverV2> create(Session session) {
        this.session = session

        final result = new ArrayList<TraceObserverV2>(5)
        createAnsiLogObserver(result)
        createGraphObserver(result)
        createReportObserver(result)
        createTimelineObserver(result)
        createTraceFileObserver(result)
        return result
    }

    protected void createAnsiLogObserver(Collection<TraceObserverV2> result) {
        if( session.ansiLog ) {
            def observer = new AnsiLogObserver()
            session.ansiLogObserver = observer
            result << observer
        }
    }

    protected void createReportObserver(Collection<TraceObserverV2> result) {
        final opts = session.config.report as Map ?: Collections.emptyMap()
        final config = new ReportConfig(opts)
        if( config.enabled )
            result << new ReportObserver(config)
    }

    protected void createTimelineObserver(Collection<TraceObserverV2> result) {
        final opts = session.config.timeline as Map ?: Collections.emptyMap()
        final config = new TimelineConfig(opts)
        if( config.enabled )
            result << new TimelineObserver(config)
    }

    protected void createGraphObserver(Collection<TraceObserverV2> result) {
        final opts = session.config.dag as Map ?: Collections.emptyMap()
        final config = new DagConfig(opts)
        if( config.enabled )
            result << new GraphObserver(config)
    }

    protected void createTraceFileObserver(Collection<TraceObserverV2> result) {
        final opts = session.config.trace as Map ?: Collections.emptyMap()
        final config = new TraceConfig(opts)
        if( config.enabled )
            result << new TraceFileObserver(config)
    }

}
