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

package nextflow.cli
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.SysEnv
import nextflow.cache.CacheDB
import nextflow.cache.CacheFactory
import nextflow.exception.AbortOperationException
import nextflow.util.HistoryFile

/**
 * Common cache operations shared by {@link CmdLog} and {@link CmdClean}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
trait CacheBase {

    @PackageScope
    Path basePath

    @PackageScope
    HistoryFile history

    abstract String getBut()

    abstract String getBefore()

    abstract String getAfter()

    abstract List<String> getArgs()

    abstract String getName()

    void init() {

        if( !history ) {
            history = !basePath ? HistoryFile.DEFAULT : new HistoryFile(basePath.resolve(HistoryFile.defaultFileName()))
        }

        if( !history.exists() || history.empty() )
            throw new AbortOperationException("It looks like no pipeline was executed in this folder (or execution history is empty)" + cloudCacheHint())

        if( after && before )
            throw new AbortOperationException("Options `after` and `before` cannot be used in the same command")

        if( after && but )
            throw new AbortOperationException("Options `after` and `but` cannot be used in the same command")

        if( before && but )
            throw new AbortOperationException("Options `before` and `but` cannot be used in the same command")

    }

    CacheDB openCache(HistoryFile.Record entry) {
        try {
            return CacheFactory.create(entry.sessionId, entry.runName, basePath).openForRead()
        }
        catch( AbortOperationException e ) {
            final hint = cloudCacheHint()
            if( !hint )
                throw e
            throw new AbortOperationException(e.message + hint, e)
        }
    }

    /**
     * Explain the failure when the cloud cache is enabled, since
     * this command can only read the local cache
     */
    private String cloudCacheHint() {
        return SysEnv.get('NXF_CLOUDCACHE_PATH')
            ? "\n\nNote: the `${getName()}` command does not support the cloud cache"
            : ''
    }

    List<HistoryFile.Record> listIds() {

        if( but ) {
            return history.findBut(but)
        }

        if( before ) {
            return history.findBefore(before)
        }

        else if( after ) {
            return history.findAfter(after)
        }

        // -- get the session ID from the command line if specified or retrieve from
        if( !args )
            return history.findByIdOrName('last')

        List<HistoryFile.Record> result = []
        for( String name : args ) {
            result.addAll(history.findByIdOrName(name))
        }
        return result
    }


}
