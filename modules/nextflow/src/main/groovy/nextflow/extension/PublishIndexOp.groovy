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
import nextflow.Global
import nextflow.Session
import nextflow.util.CsvWriter
/**
 * Publish an index file describing all files from a source
 * channel, including metadata.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class PublishIndexOp {

    private DataflowReadChannel source

    private Path basePath

    private Path path

    private Closure mapper

    private /* boolean | List<String> */ header = false

    private String sep = ','

    private List records = []

    private Session getSession() { Global.session as Session }

    PublishIndexOp(DataflowReadChannel source, Path basePath, String indexPath, Map opts) {
        this.source = source
        this.basePath = basePath
        this.path = basePath.resolve(indexPath)
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
        final record = mapper != null ? mapper.call(value) : value
        final normalized = normalizePaths(record)
        log.trace "Normalized record for index file: ${normalized}"
        records << normalized
    }

    protected void onComplete(nope) {
        if( records.size() == 0 )
            return
        log.trace "Saving records to index file: ${records}"
        new CsvWriter(header: header, sep: sep).apply(records, path)
    }

    protected Object normalizePaths(value) {
        if( value instanceof Collection ) {
            return value.collect { el ->
                if( el instanceof Path )
                    return normalizePath(el)
                if( el instanceof Collection<Path> )
                    return normalizePaths(el)
                return el
            }
        }

        if( value instanceof Map ) {
            return value.collectEntries { k, v ->
                if( v instanceof Path )
                    return List.of(k, normalizePath(v))
                if( v instanceof Collection<Path> )
                    return List.of(k, normalizePaths(v))
                return List.of(k, v)
            }
        }

        throw new IllegalArgumentException("Index file record must be a list or map: ${value} [${value.class.simpleName}]")
    }

    private Path normalizePath(Path path) {
        final sourceDir = getTaskDir(path)
        return basePath.resolve(sourceDir.relativize(path))
    }

    /**
     * Given a path try to infer the task directory to which the path below
     * ie. the directory starting with a workflow work dir and having at lest
     * two sub-directories eg work-dir/xx/yyyyyy/etc
     *
     * @param path
     */
    protected Path getTaskDir(Path path) {
        if( path == null )
            return null
        return getTaskDir0(path, session.workDir.resolve('tmp'))
            ?: getTaskDir0(path, session.workDir)
            ?: getTaskDir0(path, session.bucketDir)
    }

    private Path getTaskDir0(Path file, Path base) {
        if( base == null )
            return null
        if( base.fileSystem != file.fileSystem )
            return null
        final len = base.nameCount
        if( file.startsWith(base) && file.getNameCount() > len+2 )
            return base.resolve(file.subpath(len,len+2))
        return null
    }

}
