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

package nextflow.cloud.google.util

import com.google.api.gax.retrying.RetrySettings
import com.google.cloud.storage.Storage
import com.google.cloud.storage.StorageOptions
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.Global
import nextflow.Session
import nextflow.cloud.google.GoogleOpts

/**
 * Builds the GCS {@link Storage} client used by the global cache — the ranged reads, the flat
 * listings and the work-dir claim — from the configured {@link GoogleOpts}.
 *
 * <p>Exists because those three used to hand-build a client with credentials and project id only,
 * silently discarding the transport timeouts and retry policy the user configured. Note that this is
 * the CACHE's client: {@code GsPathFactory} builds its own for pipeline NIO traffic, and applies the
 * complementary half (timeouts and retry from an ADC-resolved default, no explicit credentials).
 * Unifying the two would change how all pipeline I/O authenticates, which is a decision of its own.
 *
 * <p>Takes {@code GoogleOpts} as a parameter rather than resolving the session itself, so each
 * caller keeps its own policy for a missing session.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
class GsStorageOptions {

    /** A client carrying credentials, project id, transport timeouts and the retry policy. */
    static Storage clientFor(GoogleOpts opts) {
        return optionsFor(opts).getService()
    }

    /**
     * The shared client for a set of options, built once per JVM.
     *
     * <p>Both GCS extensions of the global cache -- the object-store reader and the lock provider --
     * derive their client from the same {@link GoogleOpts}, and each memoized its own, so a run held
     * two clients (two connection pools, two credential refreshers) that could only ever be
     * identical. They cannot share a base class, since they extend different SPIs, so they share the
     * memo instead. Keyed by the options: a test that swaps the session config gets its own client
     * rather than a stale one.
     */
    @Memoized
    static Storage sharedClientFor(GoogleOpts opts) {
        return clientFor(opts)
    }

    /**
     * The Google config of the current session, or {@code null} when there is none.
     *
     * <p>Lives here rather than on each SPI implementation because both of them -- the object-store
     * reader and the lock provider -- need exactly this and cannot share a base class, extending
     * different SPIs. {@code GoogleOpts.fromSession} dereferences the session's config unguarded, so
     * this is what makes the null-tolerance its callers assume real rather than notional.
     *
     * <p>Deliberately NOT memoized here: {@code fromSession} already is, keyed on the session, and a
     * memo on a no-arg static would outlive the session and hand a stale config to the next one.
     */
    static GoogleOpts sessionOpts() {
        final session = (Session) Global.getSession()
        return session != null ? GoogleOpts.fromSession(session) : null
    }

    static protected StorageOptions optionsFor(GoogleOpts opts) {
        final builder = StorageOptions.newBuilder()
        if( opts?.getCredentials() != null )
            builder.setCredentials(opts.getCredentials())
        if( opts?.getProjectId() )
            builder.setProjectId(opts.getProjectId())
        if( opts != null ) {
            final transport = StorageOptions.getDefaultHttpTransportOptions().toBuilder()
            if( opts.httpConnectTimeout )
                transport.setConnectTimeout( (int) opts.httpConnectTimeout.toMillis() )
            if( opts.httpReadTimeout )
                transport.setReadTimeout( (int) opts.httpReadTimeout.toMillis() )
            builder.setTransportOptions(transport.build())
            final retry = opts.storageOpts?.retryPolicy
            if( retry ) {
                builder.setRetrySettings(StorageOptions.getDefaultRetrySettings().toBuilder()
                        .setMaxAttempts(retry.maxAttempts)
                        .setRetryDelayMultiplier(retry.multiplier)
                        .setTotalTimeout(org.threeten.bp.Duration.ofSeconds(retry.maxDelaySecs()))
                        .build())
            }
        }
        return builder.build()
    }

    /**
     * The project to bill, or {@code null}. Requester-pays cannot be set on the client at all —
     * {@code StorageOptions.Builder} has no such setter — so every request the cache makes has to
     * carry it as a per-request option. Without it, a requester-pays bucket rejects every cache
     * operation while the rest of Nextflow, which does pass it through the NIO layer, works.
     */
    static String userProject(GoogleOpts opts) {
        return opts?.enableRequesterPaysBuckets ? opts.getProjectId() : null
    }
}
