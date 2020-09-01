/*
 * Copyright (c) 2019, Seqera Labs.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
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
        try {
            final fqn = stack.join('.')
            if( fqn == 'workflow.manifest.gitmodules' && object instanceof List ) {
                writeCharSequence(object.join(','), buffer)
            }
            else
                super.writeObject(key,object,buffer)
        }
        finally {
            if(key) stack.remove(pos)
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
            log.warn "Tower request field `$key` exceeds expected size | offending value: `$seq`, size: ${seq.size()} (max: $max)"
            return seq.toString().substring(0,max)
        }
        return seq
    }

}
