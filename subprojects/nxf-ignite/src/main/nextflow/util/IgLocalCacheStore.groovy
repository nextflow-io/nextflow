/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
