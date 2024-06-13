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

package nextflow.script

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.MixOp
import nextflow.extension.PublishOp
import nextflow.file.FileHelper
/**
 * Implements the DSL for publishing workflow outputs
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class OutputDsl {

    private Map<String,Map> publishConfigs = [:]

    private Path directory

    private Map defaults = [:]

    private volatile List<PublishOp> ops = []

    void directory(String directory) {
        if( this.directory )
            throw new ScriptRuntimeException("Publish directory cannot be defined more than once in the workflow publish definition")
        this.directory = FileHelper.toCanonicalPath(directory)
    }

    void contentType(String value) {
        setDefault('contentType', value)
    }

    void contentType(boolean value) {
        setDefault('contentType', value)
    }

    void ignoreErrors(boolean value) {
        setDefault('ignoreErrors', value)
    }

    void mode(String value) {
        setDefault('mode', value)
    }

    void overwrite(boolean value) {
        setDefault('overwrite', value)
    }

    void overwrite(String value) {
        setDefault('overwrite', value)
    }

    void storageClass(String value) {
        setDefault('storageClass', value)
    }

    void tags(Map value) {
        setDefault('tags', value)
    }

    void enabled( boolean value ) {
        setDefault('enabled', value)
    }

    private void setDefault(String name, Object value) {
        if( defaults.containsKey(name) )
            throw new ScriptRuntimeException("Default `${name}` option cannot be defined more than once in the workflow publish definition")
        defaults[name] = value
    }

    void target(String name, Closure closure) {
        if( publishConfigs.containsKey(name) )
            throw new ScriptRuntimeException("Target '${name}' is defined more than once in the workflow publish definition")

        final dsl = new TargetDsl()
        final cl = (Closure)closure.clone()
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.setDelegate(dsl)
        cl.call()

        publishConfigs[name] = dsl.getOptions()
    }

    void build(Map<DataflowWriteChannel,String> targets) {
        // construct mapping of target name -> source channels
        final Map<String,List<DataflowWriteChannel>> publishSources = [:]
        for( final source : targets.keySet() ) {
            final name = targets[source]
            if( !name )
                continue
            if( name !in publishSources )
                publishSources[name] = []
            publishSources[name] << source
        }

        // create publish op (and optional index op) for each target
        for( final name : publishSources.keySet() ) {
            final sources = publishSources[name]
            final mixed = sources.size() > 1
                ? new MixOp(sources.collect( ch -> CH.getReadChannel(ch) )).apply()
                : sources.first()
            final opts = publishOptions(name, publishConfigs[name] ?: [:])

            ops << new PublishOp(CH.getReadChannel(mixed), opts).apply()
        }
    }

    private Map publishOptions(String name, Map overrides) {
        if( !directory )
            directory = FileHelper.toCanonicalPath('.')

        final opts = defaults + overrides
        if( opts.containsKey('ignoreErrors') )
            opts.failOnError = !opts.remove('ignoreErrors')
        if( !opts.containsKey('overwrite') )
            opts.overwrite = 'standard'

        final path = opts.path as String ?: name
        if( path.startsWith('/') )
            throw new ScriptRuntimeException("Invalid publish target path '${path}' -- it should be a relative path")
        opts.path = directory.resolve(path)

        if( opts.index && !(opts.index as Map).path )
            throw new ScriptRuntimeException("Index file definition for publish target '${name}' is missing `path` option")

        return opts
    }

    boolean getComplete() {
        for( final op : ops )
            if( !op.complete )
                return false
        return true
    }

    static class TargetDsl {

        private Map opts = [:]

        void contentType(String value) {
            setOption('contentType', value)
        }

        void contentType(boolean value) {
            setOption('contentType', value)
        }

        void ignoreErrors(boolean value) {
            setOption('ignoreErrors', value)
        }

        void index(Closure closure) {
            final dsl = new IndexDsl()
            final cl = (Closure)closure.clone()
            cl.setResolveStrategy(Closure.DELEGATE_FIRST)
            cl.setDelegate(dsl)
            cl.call()
            setOption('index', dsl.getOptions())
        }

        void mode(String value) {
            setOption('mode', value)
        }

        void overwrite(boolean value) {
            setOption('overwrite', value)
        }

        void overwrite(String value) {
            setOption('overwrite', value)
        }

        void path(String value) {
            setOption('path', value)
        }

        void storageClass(String value) {
            setOption('storageClass', value)
        }

        void tags(Map value) {
            setOption('tags', value)
        }

        void enabled( boolean value ) {
            setOption('enabled', value)
        }

        private void setOption(String name, Object value) {
            if( opts.containsKey(name) )
                throw new ScriptRuntimeException("Publish option `${name}` cannot be defined more than once for a given target")
            opts[name] = value
        }

        Map getOptions() {
            return opts
        }

    }

    static class IndexDsl {

        private Map opts = [:]

        void header(boolean value) {
            setOption('header', value)
        }

        void header(List<String> value) {
            setOption('header', value)
        }

        void mapper(Closure value) {
            setOption('mapper', value)
        }

        void path(String value) {
            setOption('path', value)
        }

        void sep(String value) {
            setOption('sep', value)
        }

        private void setOption(String name, Object value) {
            if( opts.containsKey(name) )
                throw new ScriptRuntimeException("Index option `${name}` cannot be defined more than once for a given index definition")
            opts[name] = value
        }

        Map getOptions() {
            return opts
        }

    }

}
