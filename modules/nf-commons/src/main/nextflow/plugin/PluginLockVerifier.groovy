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

package nextflow.plugin

import java.nio.file.Files
import java.nio.file.Path
import java.security.MessageDigest

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.exception.AbortOperationException

/**
 * Verifies plugin archives against the {@code plugins.lock} file.
 *
 * The feature is opt-in <em>by the presence of the lock file</em>: it is dormant (a no-op) when no
 * {@code plugins.lock} exists. When the file is present, the sha512 of the retained plugin archive
 * is re-computed locally (no network) and compared to the committed lock entry.
 *
 * Following the {@code go.sum} / {@code package-lock.json} model, the lock is populated
 * automatically: the first time a coordinate is downloaded and it is missing from the lock, its
 * archive checksum is appended (trust-on-first-use). An existing entry is never rewritten silently
 * — a mismatch against a committed entry is a verification failure.
 *
 * The behaviour on a checksum mismatch is gated by the {@code NXF_PLUGINS_LOCK_MODE} environment
 * variable:
 * <ul>
 *   <li>{@code strict} - abort with an {@link AbortOperationException}</li>
 *   <li>{@code warn} (default) - log a warning once per coordinate and proceed</li>
 *   <li>{@code off} - skip verification silently</li>
 * </ul>
 *
 * A locked plugin whose archive is not available locally (e.g. a cache extracted before this
 * feature existed) is never aborted: it cannot be verified offline and its archive is never
 * re-downloaded just to verify it. Integrity of the extracted code that actually runs is the
 * responsibility of the plugin directory guard, not of this archive check.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class PluginLockVerifier {

    static enum Mode { STRICT, WARN, OFF }

    private final Path lockPath

    private final PluginLockFile lock

    private final boolean enabled

    private final Set<String> notified = Collections.synchronizedSet(new HashSet<String>())

    /**
     * @param lockPath The {@code plugins.lock} path; the feature is enabled only when this file
     *        exists. A {@code null} or missing path leaves the verifier dormant.
     */
    PluginLockVerifier(Path lockPath) {
        this.lockPath = lockPath
        this.enabled = lockPath != null && Files.exists(lockPath)
        this.lock = PluginLockFile.read(lockPath)
    }

    /**
     * @return {@code true} when a lock file is present ie. verification/pinning is enabled
     */
    boolean isEnabled() {
        return enabled
    }

    /**
     * Resolve the verification mode from the {@code NXF_PLUGINS_LOCK_MODE} environment variable.
     *
     * @return The resolved {@link Mode}; {@link Mode#WARN} when unset or unrecognised
     */
    static Mode getMode() {
        final value = SysEnv.get('NXF_PLUGINS_LOCK_MODE')
        if( !value )
            return Mode.WARN
        switch( value.toLowerCase() ) {
            case 'strict': return Mode.STRICT
            case 'warn': return Mode.WARN
            case 'off': return Mode.OFF
            default:
                log.warn "Invalid NXF_PLUGINS_LOCK_MODE value: '$value' - using default 'warn'"
                return Mode.WARN
        }
    }

    /**
     * Verify - and, on first sight, pin - a plugin archive against the lock.
     *
     * @param fqid The plugin fully-qualified id ie. {@code id@version}
     * @param zip The retained plugin archive; may be {@code null} or missing
     */
    void verify(String fqid, Path zip) {
        if( !enabled )
            return
        final entry = lock.getEntry(fqid)
        final present = zip != null && Files.exists(zip)

        // coordinate not yet locked: pin it on first download (trust-on-first-use), never fail
        if( entry == null ) {
            if( present )
                pin(fqid, sha512(zip))
            else
                log.debug "Plugin '$fqid' is not in the plugins lock file and its archive is not available to pin"
            return
        }

        // locked, but the archive is not available (e.g. a cache created before this feature):
        // it cannot be verified offline - never abort and never re-download to verify
        if( !present ) {
            if( getMode() != Mode.OFF && notified.add(fqid) )
                log.warn "Cannot verify plugin '$fqid' against the plugins lock file - its archive is not available in the cache"
            return
        }

        // verify the retained archive against the committed checksum
        final actual = sha512(zip)
        if( actual == entry.sha512 )
            return

        final reason = "Plugin '$fqid' checksum does not match the plugins lock file\n- expected: ${entry.sha512}\n- actual  : ${actual}"
        switch( getMode() ) {
            case Mode.OFF:
                log.debug "Plugins lock verification failed (ignored, mode=off) - $reason"
                break
            case Mode.STRICT:
                throw new AbortOperationException("$reason\n- delete the entry from the plugins lock file and re-run to re-pin, or restore the expected plugin archive")
            case Mode.WARN:
                // warn only once per coordinate to avoid log spam
                if( notified.add(fqid) )
                    log.warn "$reason\n- delete the entry from the plugins lock file and re-run to re-pin"
                break
        }
    }

    /**
     * Append a new entry to the lock and persist it (trust-on-first-use). Existing entries are
     * never overwritten by this path — {@link #verify} routes an already-locked coordinate through
     * the checksum comparison instead.
     */
    private synchronized void pin(String fqid, String sha512) {
        lock.addEntry(fqid, new PluginLockFile.Entry(sha512))
        if( lockPath != null )
            lock.write(lockPath)
        log.info "Added plugin '$fqid' to the plugins lock file"
    }

    /**
     * Compute the sha512 hex digest of the given file.
     *
     * @param file The file to hash
     * @return The lowercase hex-encoded sha512 digest (128 chars)
     */
    static String sha512(Path file) {
        final md = MessageDigest.getInstance('SHA-512')
        try (InputStream is = Files.newInputStream(file)) {
            final buffer = new byte[8192]
            int read
            while( (read = is.read(buffer)) != -1 )
                md.update(buffer, 0, read)
        }
        return HexFormat.of().formatHex(md.digest())
    }
}
