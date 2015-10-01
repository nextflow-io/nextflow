/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

package nextflow.util

import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.extension.FilesEx
import org.gridgain.grid.GridException
import org.gridgain.grid.cache.GridCacheTx
import org.gridgain.grid.cache.store.GridCacheStore
import org.gridgain.grid.lang.GridBiInClosure
import org.jetbrains.annotations.Nullable
/**
 *  Save the cache entries in the local file system
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GgLocalCacheStore<K,V> implements GridCacheStore<K, V> {

    private Path cachePath

    public GgLocalCacheStore() {
        cachePath = Paths.get('cache')
        FilesEx.mkdirs(cachePath)
    }

    /* USED FOR TESTING PURPOSE */
    protected GgLocalCacheStore(Session session) {
    }


    @Override
    V load(@Nullable GridCacheTx tx, K key) throws GridException {
        log.trace "Local cache load > key: $key"
        return null
    }

    @Override
    void loadCache(GridBiInClosure<K, V> clo, @Nullable Object... args) throws GridException {

    }

    @Override
    void loadAll(@Nullable GridCacheTx tx, Collection<? extends K> keys, GridBiInClosure<K, V> c) throws GridException {

    }

    @Override
    void put(@Nullable GridCacheTx tx, K key, V val) throws GridException {

    }

    @Override
    void putAll(@Nullable GridCacheTx tx, Map<? extends K, ? extends V> map) throws GridException {

    }

    @Override
    void remove(@Nullable GridCacheTx tx, K key) throws GridException {

    }

    @Override
    void removeAll(@Nullable GridCacheTx tx, Collection<? extends K> keys) throws GridException {

    }

    @Override
    void txEnd(GridCacheTx tx, boolean commit) throws GridException {

    }
}
