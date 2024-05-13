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
import nextflow.processor.PublishDir
/**
 * Publish files from a source channel.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class PublishOp {

    private DataflowReadChannel source

    private PublishDir publisher

    private volatile boolean complete

    private Session getSession() { Global.session as Session }

    PublishOp(DataflowReadChannel source, Map opts) {
        this.source = source
        this.publisher = PublishDir.create(opts)
    }

    boolean getComplete() { complete }

    PublishOp apply() {
        final events = new HashMap(2)
        events.onNext = this.&onNext
        events.onComplete = this.&onComplete
        DataflowHelper.subscribeImpl(source, events)
        return this
    }

    protected void onNext(value) {
        log.trace "Publish operator received: $value"
        final result = collectFiles([:], value)
        for( final entry : result ) {
            final sourceDir = entry.key
            final files = entry.value
            publisher.apply(files, sourceDir)
        }
    }

    protected void onComplete(nope) {
        log.trace "Publish operator complete"
        this.complete = true
    }

    protected Map<Path,Set<Path>> collectFiles(Map<Path,Set<Path>> result, value) {
        if( value instanceof Path ) {
            final sourceDir = getTaskDir(value)
            if( sourceDir !in result )
                result[sourceDir] = new HashSet(10)
            result[sourceDir] << value
        }
        else if( value instanceof Collection ) {
            for( final el : value )
                collectFiles(result, el)
        }
        return result
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
