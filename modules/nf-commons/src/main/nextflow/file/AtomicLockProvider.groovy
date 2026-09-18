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

package nextflow.file

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.plugin.Plugins
import org.pf4j.ExtensionPoint
import org.pf4j.PluginManager

/**
 * Scheme-keyed SPI for an atomic create-if-absent of a lock object in cloud storage.
 *
 * Implemented per cloud by the corresponding plugin (e.g. {@code s3}, {@code az}, {@code gs})
 * and discovered via the plugin system. The create is a direct conditional PUT
 * ({@code If-None-Match: *} / {@code ifGenerationMatch=0}) — it does <b>not</b> go through
 * {@code newOutputStream}/{@code CREATE_NEW}, so it has no side effect on existing
 * {@code CREATE_NEW} callers.
 *
 * <p>Deliberately kept separate from {@link ObjectStoreReader}, which carries the object store's
 * <i>read</i> capabilities. Both are scheme-keyed plugin extensions with the same discovery, but
 * their contracts are opposites: a read capability is optional and degrades (its methods default to
 * {@code null} and {@link ObjectStoreReader#lookup} returns {@code null} for an unserved scheme,
 * because the caller always has a fallback), whereas a lock cannot degrade at all — see
 * {@link #tryCreate} and {@link #lookup} below. Merging them would have to soften one of the two.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
abstract class AtomicLockProvider implements ExtensionPoint {

    /** @return {@code true} if this provider serves the given URI scheme. */
    abstract boolean canHandle(String scheme)

    /**
     * Create the lock object, if it does not exist yet.
     *
     * A {@code false} return means one thing only: the object already existed, i.e. this caller
     * lost the race for the lock. Any condition under which the lock <b>cannot</b> be attempted --
     * a path belonging to another file system provider, a storage error -- must be reported by
     * throwing (see {@link java.nio.file.ProviderMismatchException}), never by returning
     * {@code false}: the caller bumps the key and retries on a lost race, so a silent failure
     * turns into an endless retry loop.
     *
     * <p>There is deliberately no counterpart that removes the object. The work-dir claim this SPI
     * exists for is <b>never released</b>: the marker is what makes the hash answer "in use" across
     * runs, including after a run dies, so reclaiming it means deleting the whole work directory
     * (which the cache's clean command does) rather than releasing a lock.
     *
     * @return {@code true} iff this caller created the object.
     */
    abstract boolean tryCreate(Path lockPath)

    private static volatile List<AtomicLockProvider> providers

    /**
     * The manager the memo above was resolved against. Comparing it makes a stop/init cycle
     * invalidate the memo by itself: without it, an empty list cached while the system was up would
     * survive a restart and hide extensions the new manager does have.
     */
    private static volatile PluginManager resolvedWith

    /** Test/registration seam: inject the set of providers (bypasses plugin discovery). */
    static void setProviders(List<AtomicLockProvider> list) { providers = list; resolvedWith = Plugins.getManager() }

    static List<AtomicLockProvider> getProviders() {
        final cached = providers
        if( cached != null && resolvedWith === Plugins.getManager() )
            return cached
        List<AtomicLockProvider> result
        try {
            // priority-ordered, like every other SPI here: no implementation declares @Priority, so
            // the set is unchanged and only the iteration order becomes deterministic
            result = new ArrayList<AtomicLockProvider>(Plugins.getPriorityExtensions(AtomicLockProvider))
        }
        catch( Throwable e ) {
            log.debug "Unable to load AtomicLockProvider extensions -- Cause: ${e.message}"
            result = Collections.<AtomicLockProvider>emptyList()
        }
        // Memoize a NON-EMPTY discovery only. An empty one is ambiguous -- "no plugin implements
        // this" and "the plugins are not up yet" look identical -- and caching the second would
        // poison the JVM for the rest of the run: every claim would then throw for a scheme whose
        // provider does exist. Re-querying instead costs nothing where it matters: the providers ship
        // in the cloud plugins, and only a cache addressed by a cloud URI claims a work dir, so the
        // answer is non-empty from the first task and latches there. A FAILED lookup (the catch
        // above) is likewise never cached: that one can still resolve on a retry.
        if( result ) {
            providers = result
            resolvedWith = Plugins.getManager()
        }
        return result
    }

    /**
     * Resolve the provider registered for the given path's URI scheme.
     *
     * @return the {@link AtomicLockProvider} handling the scheme.
     * @throws IllegalStateException when no provider handles the scheme.
     */
    static AtomicLockProvider lookup(Path path) {
        final scheme = path.getFileSystem().provider().getScheme()
        final p = getProviders().find { it.canHandle(scheme) }
        if( p == null )
            throw new IllegalStateException("No AtomicLockProvider for scheme '${scheme}'")
        return p
    }

}
