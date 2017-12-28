/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
