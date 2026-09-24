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

package nextflow.cache

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.plugin.Plugins
import org.pf4j.ExtensionPoint

/**
 * Factory class that create an instance of the {@link CacheDB}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class CacheFactory implements ExtensionPoint {

    /**
     * Build the cache instance.
     *
     * <p>A factory MAY write back into the current session while resolving its cache — a
     * content-addressable cache whose work directory <b>is</b> the cache assigns
     * {@code session.workDir}, and turns {@code resumeMode} on because its hits are keyed by task
     * hash rather than by session id. This is why {@code Session.init} creates the cache before it
     * reads {@code workDir}: everything downstream of that point (the work-dir creation, the
     * observers, the {@code WorkflowMetadata} snapshot) must see the effective value, or
     * {@code workflow.workDir} reports a directory the tasks never use.
     *
     * <p>{@code SessionTest} locks that ordering, so a reorder fails a test rather than silently
     * producing a run whose reported work dir is not the one in use.
     */
    protected abstract CacheDB newInstance(UUID uniqueId, String runName, Path home=null)

    /**
     * Whether this factory can serve the current session.
     *
     * <p>Selection is by {@link nextflow.plugin.Priority} alone, so the highest-priority factory on
     * the classpath serves every session whether or not it was configured. A factory that a plugin
     * contributes for one feature therefore cannot be loaded without taking over the cache: its only
     * options are to serve a cache it was not asked for, or to abort the run. Returning {@code false}
     * is the third one — decline, and let the next factory serve.
     *
     * <p>Defaults to {@code true}, so a factory that does not override this keeps the behaviour it
     * has always had.
     *
     * @return {@code true} if this factory should be used, {@code false} to defer to the next one.
     */
    protected boolean isEnabled() { true }

    static CacheDB create(UUID uniqueId, String runName, Path home=null) {
        final factory = select(Plugins.getPriorityExtensions(CacheFactory))
        log.debug "Using Nextflow cache factory: ${factory.getClass().getName()}"
        return factory.newInstance(uniqueId, runName, home)
    }

    /**
     * The first factory that claims the session, in priority order.
     */
    @PackageScope
    static CacheFactory select(List<CacheFactory> all) {
        if( !all )
            throw new IllegalStateException("Unable to find Nextflow cache factory")
        final factory = all.find { it.isEnabled() }
        if( !factory )
            throw new IllegalStateException("Unable to find an enabled Nextflow cache factory -- tried: ${all.collect { it.getClass().getName() }.join(', ')}")
        return factory
    }

}
