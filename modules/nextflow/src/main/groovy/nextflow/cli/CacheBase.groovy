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

package nextflow.cli
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.CacheDB
import nextflow.exception.AbortOperationException
import nextflow.util.HistoryFile

import static nextflow.util.HistoryFile.Record

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

    void init() {

        if( !history ) {
            history = !basePath ? HistoryFile.DEFAULT : new HistoryFile(basePath.resolve(HistoryFile.FILE_NAME))
        }

        if( !history.exists() || history.empty() )
            throw new AbortOperationException("It looks no pipeline was executed in this folder (or execution history is empty)")

        if( after && before )
            throw new AbortOperationException("Options `after` and `before` cannot be used in the same command")

        if( after && but )
            throw new AbortOperationException("Options `after` and `but` cannot be used in the same command")

        if( before && but )
            throw new AbortOperationException("Options `before` and `but` cannot be used in the same command")

    }

    CacheDB cacheFor(Record entry) {
        new CacheDB(entry,basePath)
    }

    List<Record> listIds() {

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

        def result = []
        for( String name : args ) {
            result.addAll(history.findByIdOrName(name))
        }
        return result
    }


}
