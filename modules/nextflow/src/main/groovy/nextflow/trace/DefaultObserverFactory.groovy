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

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
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
        createAgentLogObserver(result)
        createGraphObserver(result)
        createReportObserver(result)
        createTimelineObserver(result)
        createTraceFileObserver(result)
        return result
    }

    protected void createAnsiLogObserver(Collection<TraceObserverV2> result) {
        if( session.ansiLog ) {
            def observer = new AnsiLogObserver()
            session.logObserver = observer
            result << observer
        }
    }

    protected void createAgentLogObserver(Collection<TraceObserverV2> result) {
        if( session.agentLog ) {
            def observer = new AgentLogObserver()
            observer.setStatsObserver(session.statsObserver)
            session.logObserver = observer
            result << observer
        }
    }

    protected void createReportObserver(Collection<TraceObserverV2> result) {
        final opts = session.config.report as Map ?: Collections.emptyMap()
        final config = new ReportConfig(opts)
        if( config.enabled ) {
            final path = resolve('report', config.directory, config.file)
            result << new ReportObserver(config, path)
        }
    }

    protected void createTimelineObserver(Collection<TraceObserverV2> result) {
        final opts = session.config.timeline as Map ?: Collections.emptyMap()
        final config = new TimelineConfig(opts)
        if( config.enabled ) {
            final path = resolve('timeline', config.directory, config.file)
            result << new TimelineObserver(config, path)
        }
    }

    protected void createGraphObserver(Collection<TraceObserverV2> result) {
        final opts = session.config.dag as Map ?: Collections.emptyMap()
        final config = new DagConfig(opts)
        if( config.enabled ) {
            final path = resolve('dag', config.directory, config.file)
            result << new GraphObserver(config, path)
        }
    }

    protected void createTraceFileObserver(Collection<TraceObserverV2> result) {
        final opts = session.config.trace as Map ?: Collections.emptyMap()
        final config = new TraceConfig(opts)
        if( config.enabled ) {
            final path = resolve('trace', config.directory, config.file)
            result << new TraceFileObserver(config, path)
        }
    }

    /**
     * Resolve the path of a report file. When `<scope>.directory` is given, the
     * file is resolved against it, which in turn is resolved against the workflow
     * output directory. Otherwise the file is resolved against the launch directory
     * (the legacy behavior).
     */
    protected Path resolve(String scope, String directory, String file) {
        final path = FileHelper.asPath(file)
        if( !directory )
            return path
        if( path.isAbsolute() )
            throw new AbortOperationException("Config option `${scope}.directory` cannot be used when `${scope}.file` is an absolute path -- offending value: ${file}")
        final dir = FileHelper.asPath(directory)
        if( dir.isAbsolute() )
            throw new AbortOperationException("Config option `${scope}.directory` must be a relative path -- offending value: ${directory}")
        return session.outputDir.resolve(directory).resolve(file).normalize()
    }

}
