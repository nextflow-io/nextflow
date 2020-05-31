/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class PublishOp {

    private DataflowReadChannel source

    private Map opts

    private PublishDir publisher

    private Path sourceDir

    private volatile boolean complete

    private Session getSession() { Global.session as Session }

    PublishOp(DataflowReadChannel source, Map opts) {
        this.source = source
        this.opts = opts ? new LinkedHashMap(opts) : Collections.emptyMap()

        // adapt `to`  option
        if( this.opts.containsKey('to') ) {
            this.opts.path = this.opts.to
            this.opts.remove('to')
        }

        this.publisher = PublishDir.create(this.opts)
    }

    protected boolean getComplete() { complete }

    PublishOp apply() {
        final events = new HashMap(2)
        events.onNext = this.&publish0
        events.onComplete = this.&done0
        DataflowHelper.subscribeImpl(source, events)
        return this
    }

    protected void publish0(entry) {
        log.debug "Publish operator got: $entry"
        sourceDir = null
        // use an set to avoid duplicates
        final result = new HashSet(10)
        collectFiles(entry, result)
        publisher.apply(result, sourceDir)
    }

    protected void done0(nope) {
        log.debug "Publish operator complete"
        this.complete = true
    }

    protected void collectFiles(entry, Collection<Path> result) {
        if( entry instanceof Path ) {
            result.add(entry)
            if( sourceDir == null )
                sourceDir = getTaskDir(entry)
        }
        else if( entry instanceof List ) {
            for( def x : entry ) {
                collectFiles(x, result)
            }
        }
    }

    /**
     * Given a path try to infer the task directory to which the path below
     * ie. the directory starting with a workflow work dir and having at lest
     * two sub-directories eg work-dir/xx/yyyyyy/etc
     *
     * @param path
     * @return
     */
    protected Path getTaskDir(Path path) {
        if( path==null )
            return null
        def result = getTaskDir0(path, session.workDir)
        if( result == null )
            result = getTaskDir0(path, session.bucketDir)
        return result
    }

    private Path getTaskDir0(Path file, Path base) {
        if( base==null )
            return null
        if( base.fileSystem != file.fileSystem )
            return null
        final len = base.nameCount
        if( file.startsWith(base) && file.getNameCount()>len+2 )
            return base.resolve(file.subpath(len,len+2))
        return null
    }

}
