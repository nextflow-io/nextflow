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

package nextflow.extension

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.util.CsvWriter
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class PublishIndexOp {

    private DataflowReadChannel source

    private Path path

    private Closure mapper

    private /* boolean | List<String> */ header = true

    private String sep = ','

    private List records = []

    PublishIndexOp(DataflowReadChannel source, Path path, Map opts) {
        this.source = source
        this.path = path
        if( opts.mapper )
            this.mapper = opts.mapper as Closure
        if( opts.header != null )
            this.header = opts.header
        if( opts.sep )
            this.sep = opts.sep as String
    }

    void apply() {
        final events = new HashMap(2)
        events.onNext = this.&onNext
        events.onComplete = this.&onComplete
        DataflowHelper.subscribeImpl(source, events)
    }

    protected void onNext(value) {
        final normalized = mapper != null ? mapper.call(value) : value
        log.trace "Normalized record for index file: ${normalized}"
        records << normalized
    }

    protected void onComplete(nope) {
        log.trace "Saving records to index file: ${records}"
        new CsvWriter(header: header, sep: sep).apply(records, path)
    }

}
