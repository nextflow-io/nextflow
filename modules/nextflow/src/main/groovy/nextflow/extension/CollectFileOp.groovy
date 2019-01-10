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
import static nextflow.util.CheckHelper.checkParams

import java.nio.file.Path

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Channel
import nextflow.Global
import nextflow.file.FileCollector
import nextflow.file.FileHelper
import nextflow.file.SimpleFileCollector
import nextflow.file.SortFileCollector
import nextflow.util.CacheHelper

import static nextflow.util.CacheHelper.HashMode

/**
 * Implements the body of {@link DataflowExtensions#collectFile(groovyx.gpars.dataflow.DataflowReadChannel)} operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class CollectFileOp {

    static final Map COLLECT_FILE_PARAMS = [
            sort: [Boolean,'none','true','natural','index','hash','deep',Closure,Comparator],
            seed: Object,
            name: [Path, Object],
            storeDir: [Path,File,CharSequence],
            tempDir: [Path,File,CharSequence],
            newLine: Boolean,
            sliceMaxSize: Integer,
            sliceMaxItems: Integer,
            deleteTempFilesOnClose: Boolean,
            cache: [Boolean, String],
            skip: Integer,
            keepHeader: Boolean
    ]

    private final Map params

    private DataflowQueue result

    private DataflowReadChannel channel

    private FileCollector collector

    private final Closure closure

    private Path storeDir

    private String fileName

    CollectFileOp( final DataflowReadChannel channel, Map params, final Closure closure = null ) {

        checkParams('collectFile', params, COLLECT_FILE_PARAMS)
        this.params = params
        this.channel = channel
        this.closure = closure
        this.result = new DataflowQueue()

        createFileCollector()
        defineStoreDirAndFileName()
        defineHashingParams()

        Global.onShutdown {
            // make sure to delete the collector on termination
            collector.safeClose()
        }

    }

    protected FileCollector getCollector() {
        collector
    }

    protected defineHashingParams() {

        // caching params
        collector.resumable = Global.session.resumeMode
        collector.cacheable = Global.session.cacheable && ( params?.cache?.toString() != 'false' )
        collector.hashMode = HashMode.of(params?.cache) ?: HashMode.of(Global.session.config?.process?.cache) ?: HashMode.STANDARD
        collector.hashKeys = [
                Global.session.uniqueId,
                params?.storeDir,
                params?.seed,
                params?.newLine
        ]

    }


    /*
     * If a file of an absolute path is specified, the parent
     * path is used as 'storeDir'
     */
    protected defineStoreDirAndFileName() {

        if( params?.name ) {
            if( params.name instanceof Path || params.name.toString().contains('/') ) {
                def _path = params.name as Path
                fileName = _path.name
                storeDir = _path.parent
            }
            else
                fileName = params.name
        }

        /*
         * check if a 'storeDir' is provided otherwise fallback to a temp
         * folder in the session working directory
         */
        if( params?.storeDir )
            storeDir = params?.storeDir as Path

        if( storeDir )
            storeDir.createDirIfNotExists()
        else
            storeDir = FileHelper.createTempFolder(Global.session.workDir)
    }

    /*
     * each time a value is received, invoke the closure and
     * append its result value to a file
     */
    protected processItem( item ) {
        def value = closure ? closure.call(item) : item

        // when the value is a list, the first item hold the grouping key
        // all the others values are appended
        if( value instanceof List && value.size()>1 ) {
            for( int i=1; i<value.size(); i++ ) {
                collector.add(value[0] as String, value[i])
            }
        }

        // same as above
        else if( value instanceof Object[] && value.size()>1 ) {
            for( int i=1; i<value.size(); i++ ) {
                collector.add(value[0] as String, value[i])
            }
        }

        // Path object
        else if( value instanceof Path ) {
            if( fileName )
                collector.add(fileName, value)
            else
                collector.add(value.getName(), value)
        }

        // as above
        else if( value instanceof File ) {
            if( fileName )
                collector.add(fileName, value)
            else
                collector.add(value.getName(), value)
        }

        else if( value != null ) {
            if( !fileName ) fileName = 'collect-file.data'
            collector.add(fileName, value)
        }
    }


    /*
     * Emits the files when all values have been collected
     *
     * @params obj: NOT USED. It needs to be declared because this method is invoked as a closure
     */
    protected emitItems( obj ) {
        // emit collected files to 'result' channel
        collector.saveTo(storeDir).each {
            result.bind(it)
        }
        // close the channel
        result.bind(Channel.STOP)
        // close the collector
        collector.safeClose()
    }


    protected FileCollector createFileCollector() {

        // when sorting is not required 'none' use unsorted collector
        if( params?.sort == 'none' || params?.sort == false ) {
            collector = new SimpleFileCollector()
        }
        else {
            collector = new SortFileCollector()
            switch(params?.sort) {
                case true:
                case 'true':
                case 'natural':
                    collector.sort = { it -> it }
                    break

                case 'index':
                    collector.sort = null
                    break

                case null:
                case 'hash':
                    collector.sort = { CacheHelper.hasher(it).hash().asLong() }
                    break

                case 'deep':
                    collector.sort = { CacheHelper.hasher(it, CacheHelper.HashMode.DEEP).hash().asLong() }
                    break

                case Closure:
                case Comparator:
                    collector.sort = params.sort;
                    break

                default:
                    throw new IllegalArgumentException("Not a valid collectFile `sort` parameter: ${params.sort}")
            }

            if( params?.sliceMaxSize )
                collector.sliceMaxSize = params.sliceMaxSize

            if( params?.sliceMaxItems )
                collector.sliceMaxItems = params.sliceMaxItems
        }

        // set other params
        collector.tempDir = params?.tempDir as Path
        collector.newLine = params?.newLine as Boolean
        collector.seed = params?.seed
        if( params?.deleteTempFilesOnClose != null )
            collector.deleteTempFilesOnClose = params.deleteTempFilesOnClose as boolean
        if( params?.skip )
            collector.skipLines = params?.skip
        if( params?.keepHeader != null )
            collector.keepHeader = params.keepHeader as boolean
        if( collector.keepHeader ) {
            // validate `seed` parameter
            if( collector.seed != null )
                throw new IllegalArgumentException("Parameter `keepHeader` and `seed` conflict -- check operator `collectFile`")
            // validate `skip` parameter
            if( params?.skip == null )
                collector.skipLines = 1
            else if( collector.skipLines < 1 )
                throw new IllegalArgumentException("Parameter `skip` must be greater than zero when `keepHeader` is specified -- check operator `collectFile`")
        }
        return collector
    }


    DataflowQueue apply() {
        DataflowHelper.subscribeImpl( channel, [onNext: this.&processItem, onComplete: this.&emitItems] )
        return result
    }
}
