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

package nextflow.util
import javax.cache.Cache
import javax.cache.integration.CacheLoaderException
import javax.cache.integration.CacheWriterException
import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.extension.FilesEx
import org.apache.ignite.cache.store.CacheStore
import org.apache.ignite.lang.IgniteBiInClosure
import org.jetbrains.annotations.Nullable

/**
 *  Save the cache entries in the local file system
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class IgLocalCacheStore<K,V> implements CacheStore<K, V> {

    private Path cachePath

    public IgLocalCacheStore() {
        cachePath = Paths.get('cache')
        FilesEx.mkdirs(cachePath)
    }

    /* USED FOR TESTING PURPOSE */
    protected IgLocalCacheStore(Session session) {
    }


    @Override
    V load(K key) throws CacheLoaderException {
        log.trace "Local cache load > key: $key"
        return null
    }

    @Override
    void loadCache(IgniteBiInClosure<K, V> clo, @Nullable Object... args) throws CacheLoaderException {

    }

    @Override
    void sessionEnd(boolean commit) throws CacheWriterException {

    }

    @Override
    Map<K, V>  loadAll(Iterable<? extends K> keys) throws CacheLoaderException {

    }

    @Override
    void write(Cache.Entry<? extends K, ? extends V> entry) throws CacheWriterException {

    }

    @Override
    void writeAll(Collection<Cache.Entry<? extends K, ? extends V>> entries) throws CacheWriterException {

    }

    @Override
    void delete(Object key) throws CacheWriterException {

    }

    @Override
    void deleteAll(Collection<?> keys) throws CacheWriterException {

    }
}
