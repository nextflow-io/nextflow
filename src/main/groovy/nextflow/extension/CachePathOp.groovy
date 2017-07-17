/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

package nextflow.extension
import static nextflow.util.CheckHelper.checkParams

import java.nio.file.Path

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Channel
import nextflow.Global
import nextflow.file.FileHelper
import nextflow.util.CacheHelper

/**
 * Implements the body of {@link DataflowExtensions#cachePath(groovyx.gpars.dataflow.DataflowReadChannel)} operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
@Slf4j
class CachePathOp {

    static final Map CACHE_PATH_PARAMS = [
            storeDir: [Path,File,CharSequence],
            copyLocal: Boolean,
            deep: Boolean,
    ]

    private final Map params

    private DataflowQueue result

    private DataflowReadChannel channel

    private Path storeDir

    private boolean linkLocalPaths = true

    private CacheHelper.HashMode cachingMode


    CachePathOp( final DataflowReadChannel channel, Map params) {

        checkParams('cachePath', params, CACHE_PATH_PARAMS)
        this.params = params
        this.channel = channel
        this.result = new DataflowQueue()

        defineLinkLocally()
        defineStoreDir()
        defineDeepCaching()
    }

    protected defineLinkLocally() {
        if (params?.copyLocal) {
            linkLocalPaths = !(params?.copyLocal as boolean)
        }
    }

    protected defineDeepCaching() {
        def deep = true
        if (params?.deep) {
           deep = params?.deep as boolean
        }

        if (deep) {
            cachingMode = CacheHelper.HashMode.DEEP
        } else {
            cachingMode = CacheHelper.HashMode.STANDARD
        }
    }

    protected defineStoreDir() {
        // check if a 'storeDir' is provided otherwise fallback to a temp
        // folder in the session working directory
        if( params?.storeDir )
            storeDir = params?.storeDir as Path

        if( storeDir )
            storeDir.createDirIfNotExists()
        else
            storeDir = FileHelper.createTempFolder(Global.session.workDir)
    }

    /*
     * Each time a value is received, check that it's a Path or File,
     * cache as directed, and then return the cached path.
     */
    protected processPath( inPath ) {
        if( inPath instanceof Path || inPath instanceof File ) {
            result.bind( getCachedPath( inPath ) )
        } else {
            throw new IllegalArgumentException("Can only cache files or directories!")
        }
    }

    /*
     * Main caching logic here.
     * If the input path doesn't exist in the cache, always cache it.
     * If the input path is on the local file system
     * then cache and we're symlinking local paths, then copy input into cache.
     * If the input path is remote, always copy into cache.
     * If the cache target already exists and is not a symlink, then hash
     * the files and only re-copy them if the hashes don't match.
     */
    protected getCachedPath( inPath ) {
        def cachedPath = storeDir + "/" + inPath

        if ( !cachedPath.exists() ) {

            if ( isLocal(inPath) && linkLocalPaths ) {
                log.info("Symlinking ${inPath} to ${cachedPath}")
                createSymlink( inPath, cachedPath )
            } else {
                log.info("Caching ${inPath} to ${cachedPath}")
                copy(inPath, cachedPath)
            }

        } else if ( !cachedPath.isLink() ) {

            def inHash = CacheHelper.hasher(inPath, cachingMode).hash().asLong()
            def cacheHash = CacheHelper.hasher(cachedPath, cachingMode).hash().asLong()

            if (inHash != cacheHash) {
                log.info("Re-Caching ${inPath} to ${cachedPath}")
                copy(inPath, cachedPath)
            }
        }

        return cachedPath
    }

    protected createSymlink(inPath, cachedPath) {
        cachedPath.getParent().mkdirs()
        inPath.mklink( cachedPath )
    }

    protected copy(inPath, cachedPath) {
        if ( cachedPath.exists() && cachedPath.isDirectory() ) {
            FileHelper.deletePath(cachedPath)
        }
        inPath.copyTo( cachedPath )
    }

    protected isLocal(path) {
        def pathUri = path.toUri().toString()
        return pathUri.startsWith("file")
    }

    protected finish(obj) {
        result.bind(Channel.STOP)
    }

    DataflowQueue apply() {
        DataflowHelper.subscribeImpl( channel, [onNext: this.&processPath, onComplete: this.&finish] )
        return result
    }
}
