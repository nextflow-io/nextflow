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
 * Scheme-keyed SPI for the two <b>read</b> operations an object store supports natively but the NIO
 * {@code FileSystem} abstraction hides:
 *
 * <ul>
 *   <li>{@link #listWithMeta} — every member of a prefix in a single flat, <b>delimiter-less</b>
 *       listing, instead of one delimiter-based LIST per subdirectory via the NIO walk;</li>
 *   <li>{@link #readRange} — a byte range as a single ranged GET, instead of downloading (or
 *       streaming from the start of) the whole object.</li>
 * </ul>
 *
 * Used by the record-backed file identities: the directory identity is folded from one listing, and
 * the {@code sample} identity reads a few small windows of a large object. Both are <b>optional</b>:
 * each method returns {@code null} by default and {@link #lookup} returns {@code null} for a scheme
 * no provider handles, so a caller always has a documented fallback (the hierarchical walk, a full
 * read, or default path-based hashing).
 *
 * <p>Deliberately <b>not</b> merged with {@link AtomicLockProvider}, whose contract is the opposite
 * on every point that matters: a conditional create cannot degrade gracefully (a {@code false} that
 * did not mean "lost the race" would spin the caller's bump-and-retry forever), so its methods are
 * abstract and its {@code lookup} throws rather than returning {@code null}. The atomicity guarantee
 * lives in that type's name, and one shared {@code lookup} could not answer for both.
 *
 * <p>Declared as an abstract class rather than an interface: the static resolver below is called from
 * classes that are not {@code @CompileStatic}, and Groovy's dynamic dispatch cannot invoke a static
 * method declared on an interface.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
abstract class ObjectStoreReader implements ExtensionPoint {

    private static volatile List<ObjectStoreReader> providers

    /**
     * The manager the memo above was resolved against. Comparing it makes a stop/init cycle
     * invalidate the memo by itself: without it, an empty list cached while the system was up would
     * survive a restart and hide extensions the new manager does have.
     */
    private static volatile PluginManager resolvedWith

    /** @return {@code true} if this provider handles the given URI scheme (e.g. {@code s3}). */
    abstract boolean canHandle(String scheme)

    /**
     * Flat-list every object under {@code prefix} (a directory) as {@code (relpath-relative-to-prefix,
     * {@link ObjectMeta})} pairs, using a single paginated LIST with <b>no delimiter</b> — so the cost
     * is {@code ceil(N/1000)} calls regardless of subfolder nesting, versus one delimiter-based LIST
     * per subdirectory via the NIO walk. Returns {@code null} by default (not implemented); callers
     * then fall back to a hierarchical walk.
     */
    List<Map.Entry<String,ObjectMeta>> listWithMeta(Path prefix) { return null }

    /**
     * Read exactly {@code len} bytes starting at {@code offset} from the object at {@code path} as a
     * single ranged GET. Returns {@code null} by default (not implemented); callers then fall back to
     * a full (non-ranged) read.
     */
    byte[] readRange(Path path, long offset, int len) { return null }

    /**
     * Normalize a directory prefix so it ends with a single {@code '/'} (empty stays empty, i.e. the
     * bucket root). Shared by implementations so the LIST prefix is derived identically on every cloud.
     */
    protected static String normalizePrefix(String key) {
        return !key ? '' : (key.endsWith('/') ? key : key + '/')
    }

    /**
     * The relpath of an object {@code name} relative to the (already-normalized) {@code base} prefix,
     * or {@code null} to skip it — a "directory" placeholder object (trailing {@code '/'}) or the base
     * prefix itself. Shared so directory members are filtered/relativized identically on every cloud
     * (the identities must match across providers).
     */
    protected static String relativize(String base, String name) {
        if( name.endsWith('/') )
            return null
        final rel = base ? name.substring(base.length()) : name
        return rel ?: null
    }

    /**
     * Test/registration seam: inject the set of providers (bypasses plugin discovery). Public, not
     * package-scoped: the consumers that need to inject a fake live in other modules.
     */
    static void setProviders(List<ObjectStoreReader> it) { providers = it; resolvedWith = Plugins.getManager() }

    static List<ObjectStoreReader> getProviders() {
        final cached = providers
        if( cached != null && resolvedWith === Plugins.getManager() )
            return cached
        try {
            final result = new ArrayList<ObjectStoreReader>(Plugins.getPriorityExtensions(ObjectStoreReader))
            // Memoize a NON-EMPTY discovery only. An empty one is ambiguous -- "no plugin implements
            // this" and "the plugins are not up yet" look identical -- and caching the second would
            // poison the JVM for the rest of the run, silently disabling the identity, since `lookup`
            // degrades to null. Re-querying instead costs nothing where it matters: the readers ship
            // in the cloud plugins, and a cache addressed by a cloud URI has one loaded, so the
            // answer is non-empty from the first task and latches there. A FAILED lookup (the catch
            // below) is likewise never cached: it can still resolve later.
            if( result ) {
                providers = result
                resolvedWith = Plugins.getManager()
            }
            return result
        }
        catch( Throwable e ) {
            log.debug "Unable to load the object-store readers -- ${e.message}"
            return Collections.<ObjectStoreReader>emptyList()
        }
    }

    /** Resolve the provider for the path's scheme, or {@code null} if none handles it. */
    static ObjectStoreReader lookup(Path path) {
        final scheme = path.getFileSystem().provider().getScheme()
        return getProviders().find { it.canHandle(scheme) }
    }
}
