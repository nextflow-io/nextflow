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

import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes
import java.security.MessageDigest

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.exception.AbortOperationException

/**
 * Verifies an extracted plugin directory against the {@code plugins.lock} file.
 *
 * The feature is opt-in <em>by the presence of the lock file</em>: it is dormant (a no-op) when no
 * {@code plugins.lock} exists. When the file is present, a canonical hash of the extracted
 * {@code $id-$version/} directory - the code that Nextflow actually loads and executes - is
 * re-computed locally (no network) and compared to the committed lock entry. Because it hashes the
 * unpacked tree rather than the archive, it detects both a tampered/compromised download and a
 * poisoned cache directory (a lower-trust user editing already-extracted files on a shared cache).
 *
 * Following the {@code go.sum} / {@code package-lock.json} model the lock is populated
 * automatically: the first time a coordinate is seen and it is missing from the lock, its tree
 * hash is appended (trust-on-first-use). An existing entry is never rewritten silently — a
 * mismatch against a committed entry is a verification failure, gated by
 * {@code NXF_PLUGINS_LOCK_MODE}:
 * <ul>
 *   <li>{@code strict} - abort with an {@link AbortOperationException}</li>
 *   <li>{@code warn} (default) - log a warning once per coordinate and proceed</li>
 *   <li>{@code off} - disable the feature altogether ie. no hashing, no pinning, no reporting</li>
 * </ul>
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
     *        exists. A {@code null} or missing path leaves the verifier dormant, and so does
     *        {@code NXF_PLUGINS_LOCK_MODE=off}.
     */
    PluginLockVerifier(Path lockPath) {
        this.lockPath = lockPath
        // `off` is a complete escape hatch: the lock file is not even read, so a run can always
        // be unblocked without deleting a file committed in the pipeline repository
        this.enabled = lockPath != null && Files.exists(lockPath) && getMode() != Mode.OFF
        this.lock = enabled ? PluginLockFile.read(lockPath) : new PluginLockFile()
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
     * Verify - and, on first sight, pin - an extracted plugin directory against the lock.
     *
     * @param fqid The plugin fully-qualified id ie. {@code id@version}
     * @param pluginDir The extracted plugin directory ({@code $id-$version/})
     */
    void verify(String fqid, Path pluginDir) {
        if( !enabled )
            return
        if( pluginDir == null || !Files.isDirectory(pluginDir) ) {
            log.debug "Cannot verify plugin '$fqid' against the plugins lock file - directory not available: $pluginDir"
            return
        }

        final actual = sha512Tree(pluginDir)
        final entry = lock.getEntry(fqid)

        // coordinate not yet locked: pin it on first sight (trust-on-first-use), never fail
        if( entry == null ) {
            pin(fqid, actual)
            return
        }
        // matches the committed hash
        if( actual == entry.sha512 )
            return

        // mismatch: the extracted plugin differs from what the lock pinned
        // note: `off` is handled upfront by the `enabled` flag, so only strict/warn get here
        final reason = "Plugin '$fqid' does not match the plugins lock file\n- expected: ${entry.sha512}\n- actual  : ${actual}"
        if( getMode() == Mode.STRICT )
            throw new AbortOperationException("$reason\n- delete the entry from the plugins lock file and re-run to re-pin, or restore the expected plugin")
        // warn only once per coordinate to avoid log spam
        if( notified.add(fqid) )
            log.warn "$reason\n- delete the entry from the plugins lock file and re-run to re-pin"
    }

    /**
     * Append a new entry to the lock and persist it (trust-on-first-use). Existing entries are
     * never overwritten by this path — {@link #verify} routes an already-locked coordinate through
     * the hash comparison instead.
     *
     * Failing to persist the pin is not an integrity problem, therefore it is reported as a
     * warning instead of aborting the run eg. when the lock file lives in a read-only checkout.
     */
    private synchronized void pin(String fqid, String sha512) {
        lock.addEntry(fqid, new PluginLockFile.Entry(sha512))
        try {
            lock.write(lockPath)
            log.info "Added plugin '$fqid' to the plugins lock file"
        }
        catch( Exception e ) {
            log.warn "Unable to add plugin '$fqid' to the plugins lock file: $lockPath - cause: ${e.message ?: e}"
        }
    }

    /**
     * Compute a canonical sha512 digest over the content of a directory tree. The digest covers,
     * for every regular file (visited in sorted relative-path order for determinism), its relative
     * path and its bytes - not timestamps or permissions - so it is stable across extractions and
     * platforms while still detecting any change to the files that will be executed.
     *
     * @param dir The directory to hash
     * @return The lowercase hex-encoded sha512 digest (128 chars)
     */
    static String sha512Tree(Path dir) {
        final md = MessageDigest.getInstance('SHA-512')
        // TreeMap keyed by relative path -> deterministic, sorted iteration order
        final files = new TreeMap<String, Path>()
        Files.walkFileTree(dir, new SimpleFileVisitor<Path>() {
            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
                // only regular files are hashed; symlinks and other special entries are skipped
                // because they cannot be read as a byte stream in a portable way
                if( attrs.isRegularFile() )
                    files.put(dir.relativize(file).toString().replace('\\', '/'), file)
                return FileVisitResult.CONTINUE
            }
        })
        final buffer = new byte[8192]
        for( Map.Entry<String, Path> it : files.entrySet() ) {
            md.update(it.key.getBytes('UTF-8'))
            md.update((byte) 0)
            try (InputStream is = Files.newInputStream(it.value)) {
                int read
                while( (read = is.read(buffer)) != -1 )
                    md.update(buffer, 0, read)
            }
            md.update((byte) 0)
        }
        return HexFormat.of().formatHex(md.digest())
    }
}
