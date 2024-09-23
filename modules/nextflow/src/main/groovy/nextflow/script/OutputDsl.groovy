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
import nextflow.Global
import nextflow.Session
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

    private Session session = Global.session as Session

    private Map<String,Map> targetConfigs = [:]

    private volatile List<PublishOp> ops = []

    void target(String name, Closure closure) {
        if( targetConfigs.containsKey(name) )
            throw new ScriptRuntimeException("Target '${name}' is defined more than once in the workflow output definition")

        final dsl = new TargetDsl()
        final cl = (Closure)closure.clone()
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.setDelegate(dsl)
        cl.call()

        targetConfigs[name] = dsl.getOptions()
    }

    void build(Map<DataflowWriteChannel,String> targets) {
        final defaults = session.config.navigate('workflow.output', Collections.emptyMap()) as Map

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

        // validate target configs
        for( final name : targetConfigs.keySet() ) {
            if( name !in publishSources )
                log.warn "Publish target '${name}' was defined in the output block but not used by the workflow"
        }

        // create publish op (and optional index op) for each target
        for( final name : publishSources.keySet() ) {
            final sources = publishSources[name]
            final mixed = sources.size() > 1
                ? new MixOp(sources.collect( ch -> CH.getReadChannel(ch) )).apply()
                : sources.first()
            final overrides = targetConfigs[name] ?: Collections.emptyMap()
            final opts = publishOptions(name, defaults, overrides)

            if( opts.enabled == null || opts.enabled )
                ops << new PublishOp(CH.getReadChannel(mixed), opts).apply()
        }
    }

    private Map publishOptions(String name, Map defaults, Map overrides) {
        final opts = defaults + overrides
        if( opts.containsKey('ignoreErrors') )
            opts.failOnError = !opts.remove('ignoreErrors')
        if( !opts.containsKey('overwrite') )
            opts.overwrite = 'standard'

        final path = opts.path as String ?: name
        if( path.startsWith('/') || path.endsWith('/') )
            throw new ScriptRuntimeException("Invalid publish target path '${path}' -- it should not contain a leading or trailing slash")
        opts.path = session.outputDir.resolve(path)

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

        void enabled(boolean value) {
            setOption('enabled', value)
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

        void path(Closure value) {
            setOption('path', '.')
            setOption('pathAs', value)
        }

        void storageClass(String value) {
            setOption('storageClass', value)
        }

        void tags(Map value) {
            setOption('tags', value)
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
