/*
 * Copyright 2013-2024, Seqera Labs
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

import nextflow.Session
import nextflow.file.FileHelper

/**
 * Creates Nextflow observes object
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultObserverFactory implements TraceObserverFactoryV2 {

    private Map config

    @Override
    Collection<TraceObserverV2> create(Session session) {
        this.config = session.config

        final result = new ArrayList<TraceObserverV2>(5)
        if( session.ansiLog )
            createAnsiLogObserver(result)
        createGraphObserver(result)
        createReportObserver(result)
        createTimelineObserver(result)
        createTraceFileObserver(result)
        return result
    }

    protected void createAnsiLogObserver(Collection<TraceObserverV2> result) {
        session.ansiLogObserver = new AnsiLogObserver()
        result << session.ansiLogObserver
    }

    protected void createReportObserver(Collection<TraceObserverV2> result) {
        final isEnabled = config.navigate('report.enabled') as Boolean
        if( !isEnabled )
            return

        final fileName = config.navigate('report.file', ReportObserver.DEF_FILE_NAME) as String
        final maxTasks = config.navigate('report.maxTasks', ReportObserver.DEF_MAX_TASKS) as int
        final overwrite = config.navigate('report.overwrite') as Boolean

        final observer = new ReportObserver(FileHelper.asPath(fileName), overwrite)
        observer.maxTasks = maxTasks
        result << observer
    }

    protected void createTimelineObserver(Collection<TraceObserverV2> result) {
        final isEnabled = config.navigate('timeline.enabled') as Boolean
        if( !isEnabled )
            return

        final fileName = config.navigate('timeline.file', TimelineObserver.DEF_FILE_NAME) as String
        final overwrite = config.navigate('timeline.overwrite') as Boolean

        result << new TimelineObserver(FileHelper.asPath(fileName), overwrite)
    }

    protected void createGraphObserver(Collection<TraceObserverV2> result) {
        final isEnabled = config.navigate('dag.enabled') as Boolean
        if( !isEnabled )
            return

        final fileName = config.navigate('dag.file', GraphObserver.DEF_FILE_NAME) as String
        final overwrite = config.navigate('dag.overwrite') as Boolean

        result << new GraphObserver(FileHelper.asPath(fileName), overwrite)
    }

    protected void createTraceFileObserver(Collection<TraceObserverV2> result) {
        final isEnabled = config.navigate('trace.enabled') as Boolean
        if( !isEnabled )
            return

        final fields = config.navigate('trace.fields', '') as String
        final fileName = config.navigate('trace.file', TraceFileObserver.DEF_FILE_NAME) as String
        final overwrite = config.navigate('trace.overwrite') as Boolean
        final raw = config.navigate('trace.raw') as Boolean
        final separator = config.navigate('trace.sep', '\t') as String

        final observer = new TraceFileObserver(FileHelper.asPath(fileName), overwrite, separator)
        observer.useRawNumbers(raw)
        observer.setFieldsAndFormats(fields)
        result << observer
    }

}
