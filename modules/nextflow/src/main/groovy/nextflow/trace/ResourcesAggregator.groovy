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

import java.util.concurrent.Callable
import java.util.concurrent.Executors
import java.util.concurrent.Future

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Collect and aggregate execution metrics used by execution report
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ResourcesAggregator {

    /**
     * Holds a report summary instance for each process group
     */
    final private Map<String,ReportSummary> summaries = new LinkedHashMap<>()

    /**
     * Aggregates task record for each process in order to render the
     * final execution stats
     *
     * @param record A {@link TraceRecord} object representing a task executed
     */
    void aggregate(TraceRecord record) {
        // aggregate on the process simple name
        // therefore all nested process are kept together
        def process = record.getSimpleName()
        def summary = summaries.get(process)
        if( !summary ) {
            summaries.put(process, summary=new ReportSummary())
        }
        summary.add(record)
    }

    /**
     * Compute the workflow process summary stats
     *
     * @return A {@link Map} holding the summary stats for each process
     */
    protected Map<String,Map> computeSummaryMap() {

        final result = new LinkedHashMap<String,Map>(summaries.size())

        // summary stats can be expensive on big workflow
        // speed-up the computation using a parallel
        List<Callable<List>> tasks = []
        summaries.keySet().each  { String process ->
            final summary = summaries[process]
            summary.names.each { String series ->
                // the task execution turn a triple
                tasks << { return [ process, series, summary.compute(series)] } as Callable<List>
            }
            // initialise the result entry
            result.put(process, new HashMap(10))
        }

        // use a dedicated thread pool to avoid starvation deadlock
        // with the shared virtual thread executor (see #6833)
        final executor = Executors.newFixedThreadPool(Runtime.runtime.availableProcessors())
        try {
            // submit the parallel execution
            final allResults = executor.invokeAll(tasks)

            // compose the final result
            for( Future<List> future : allResults ) {
                final triple = future.get()
                final name = triple[0]      // the process name
                final series = triple[1]    // the series name eg. `cpu`, `time`, etc
                final summary = triple[2]   // the computed summary
                result.get(name).put(series, summary)
            }
        }
        finally {
            executor.shutdown()
        }

        return result
    }

    List<Map> computeSummaryList() {
        def map = computeSummaryMap()
        def result = new ArrayList(map.size())
        for( Map.Entry<String,Map> entry : map.entrySet() ) {
            def record = entry.value
            record.process = entry.key
            result.add( record )
        }
        return result
    }


    /**
     * @return The execution summary json
     */
    String renderSummaryJson() {
        final summary = computeSummaryList()
        return JsonOutput.toJson(summary)
    }
}
