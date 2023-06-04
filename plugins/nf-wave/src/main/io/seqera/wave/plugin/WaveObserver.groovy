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

import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import io.seqera.wave.plugin.config.ReportOpts
import nextflow.Session
import nextflow.file.FileHelper
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class WaveObserver implements TraceObserver {

    private WaveClient client

    private ConcurrentHashMap<String,String> containers = new ConcurrentHashMap<>()

    WaveObserver(Session session) {
        this.client = new WaveClient(session)
    }

    protected void apply(TaskHandler handler) {
        final process = handler.task.getProcessor().getName()
        containers.computeIfAbsent(process, (String it) -> {
            final container = handler.task.getContainer()
            return client.resolveSourceContainer(container)
        })
    }

    void onProcessComplete(TaskHandler handler, TraceRecord trace){
        apply(handler)
    }

    void onProcessCached(TaskHandler handler, TraceRecord trace){
        apply(handler)
    }

    @Override
    void onFlowComplete() {
        final result = renderContainersConfig(containers)
        // save the report file
        FileHelper
                .asPath(reportOpts().file())
                .text = result.toString()
    }

    protected String renderContainersConfig(Map<String,String> containers) {
        final result = new StringBuilder()
        for( Map.Entry<String,String> entry : containers ) {
            result.append("process { withName: '${entry.key}' { container='$entry.value' }}\n")
        }
        return result.toString()
    }

    ReportOpts reportOpts() {
        client.config().reportOpts()
    }
}
