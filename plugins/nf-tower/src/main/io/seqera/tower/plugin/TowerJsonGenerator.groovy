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
 *
 */

package io.seqera.tower.plugin

import java.nio.file.Path
import java.time.OffsetDateTime

import groovy.json.DefaultJsonGenerator
import groovy.json.JsonGenerator
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.NextflowMeta
import nextflow.trace.ProgressRecord
import nextflow.util.Duration
import org.apache.groovy.json.internal.CharBuf
/**
 * Customized json generator that chomp string values
 * longer than expected size
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerJsonGenerator extends DefaultJsonGenerator {

    List<String> stack = new ArrayList<>(10)
    Map<String,Integer> scheme

    static TowerJsonGenerator create(Map<String,Integer> scheme) {
        final opts = new JsonGenerator.Options()
                .addConverter(Path) { Path p, String key -> p.toUriString() }
                .addConverter(Duration) { Duration d, String key -> d.durationInMillis }
                .addConverter(NextflowMeta) { NextflowMeta m, String key -> m.toJsonMap() }
                .addConverter(OffsetDateTime) { it.toString() }
                .dateFormat(Const.ISO_8601_DATETIME_FORMAT).timezone("UTC")

        return new TowerJsonGenerator(opts, scheme)
    }

    protected TowerJsonGenerator(Options options, Map<String,Integer> scheme) {
        super(options)
        this.scheme = scheme
    }

    @Override
    protected Map<?, ?> getObjectProperties(Object object) {
        final result = super.getObjectProperties(object)
        if( object instanceof ProgressRecord ) {
            result.remove('hash')
            result.remove('errored')
            result.remove('completedCount')
            result.remove('totalCount')
            result.remove('taskName')
        }
        return result
    }

    @Override
    protected void writeObject(String key, Object object, CharBuf buffer) {
        final pos=stack.size()
        if(key) stack.add(pos, key)
        final fqn = stack.join('.')
        try {
            if( fqn == 'workflow.manifest.gitmodules' && object instanceof List ) {
                writeCharSequence(object.join(','), buffer)
            }
            else
                super.writeObject(key,object,buffer)
        }
        catch( Exception e ) {
            log.warn1 ("Unable to serialize key=$fqn; value=${safeString0(object)}; type=${object?.getClass()?.getName()} -- Cause: ${e.message ?: e}", causedBy: e)
        }
        finally {
            if(key) stack.remove(pos)
        }
    }

    private String safeString0(Object obj) {
        try {
            return obj !=null ? obj.toString() : null
        }
        catch( Throwable e ) {
            log.debug "SafeString error=${e.message ?: e}"
            return null
        }
    }

    @Override
    protected void writeRaw(CharSequence seq, CharBuf buffer) {
        super.writeRaw(chompString(seq),buffer)
    }

    @Override
    protected void writeCharSequence(CharSequence seq, CharBuf buffer) {
        super.writeCharSequence(chompString(seq), buffer)
    }

    /**
     * Make sure the specified argument is not longer than the expected length
     * defined in the `scheme` object
     *
     * @param seq The string object as {@link CharSequence} instance
     * @return A string object whose length does not exceed the max for the current key entry
     */
    final protected CharSequence chompString(CharSequence seq) {
        final key = stack.join('.')
        final max = scheme.get(key)
        if( seq!=null && max && seq.length()>max ) {
            final result = seq.toString().substring(0,max)
            // show only the first 100 chars in the log as a preview
            final preview = result.length()>100 ? result.substring(0,100) + '(truncated)' : result
            log.warn "Seqera Platform request field `$key` exceeds expected size | offending value: `${preview}`, size: ${seq.size()} (max: $max)"
            return result
        }
        return seq
    }

}
