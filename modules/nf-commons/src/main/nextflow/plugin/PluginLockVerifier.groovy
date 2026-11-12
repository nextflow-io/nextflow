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
 *   <li>{@code strict} - abort with an {@link AbortOperationException}, both on a hash mismatch and
 *       when a plugin reaches the loader without a corresponding lock entry (nothing is auto-pinned,
 *       so a version that resolved past the pin and is already extracted in the cache cannot slip
 *       through unverified)</li>
 *   <li>{@code warn} (default) - log a warning once per coordinate and proceed; a coordinate missing
 *       from the lock is pinned on first sight (trust-on-first-use)</li>
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
        if( entry == null )
            pinOrRefuse(fqid, actual)            // not yet locked: TOFU pin, or strict refusal
        else if( actual != entry.sha512 )
            reportMismatch(fqid, actual, entry)  // differs from the pinned hash
    }

    /**
     * Handle a coordinate not yet in the lock. Strict mode refuses it: a plugin reaching the loader
     * unpinned has bypassed verification — an unpinned coordinate can resolve past the committed
     * version to a tree already extracted in a shared cache, where the download short-circuit means
     * pf4j's own verifier never ran either. Otherwise pin it on first sight (trust-on-first-use).
     */
    private void pinOrRefuse(String fqid, String actual) {
        if( getMode() == Mode.STRICT )
            throw new AbortOperationException(
                    "Plugin '$fqid' is not present in the plugins lock file: $lockPath\n" +
                    "- pin the plugin version in the `plugins` config scope, or re-run with NXF_PLUGINS_LOCK_MODE=warn to record it")
        warnOnVersionDrift(fqid)
        pin(fqid, actual)
    }

    /**
     * Handle an extracted tree that differs from the committed hash: abort in strict mode, otherwise
     * warn once per coordinate. ({@code off} never reaches here — the verifier is disabled upfront.)
     */
    private void reportMismatch(String fqid, String actual, PluginLockFile.Entry entry) {
        final reason = "Plugin '$fqid' does not match the plugins lock file\n- expected: ${entry.sha512}\n- actual  : ${actual}"
        if( getMode() == Mode.STRICT )
            throw new AbortOperationException("$reason\n- delete the entry from the plugins lock file and re-run to re-pin, or restore the expected plugin")
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
     * Warn when the lock already pins a <em>different</em> version of the same plugin id than the
     * one about to be pinned. This usually means an unpinned plugin resolved to a newer release
     * than the one recorded in the lock, so the drift is surfaced instead of silently accumulating
     * a second entry for the same plugin.
     */
    private void warnOnVersionDrift(String fqid) {
        final p = fqid.lastIndexOf('@')
        if( p < 0 )
            return
        // include the trailing '@' so 'nf-amazon@' does not match 'nf-amazon-extra@...'
        final prefix = fqid.substring(0, p + 1)
        final others = lock.getEntries().keySet().findAll { it != fqid && it.startsWith(prefix) }
        if( others )
            log.warn "Plugin '$fqid' is being added to the plugins lock file, which already pins a different version: ${others.join(', ')}"
    }

    /**
     * Compute a canonical sha512 digest over the content of a directory tree. Visiting entries in
     * sorted relative-path order for determinism, the digest covers every regular file (its relative
     * path and its bytes) and every symbolic link (its relative path and its target) - not
     * timestamps or permissions - so it is stable across extractions and platforms while still
     * detecting any change to the files, or links, that will be loaded and executed.
     *
     * @param dir The directory to hash
     * @return The lowercase hex-encoded sha512 digest (128 chars)
     */
    static String sha512Tree(Path dir) {
        final md = MessageDigest.getInstance('SHA-512')
        for( Map.Entry<String,Path> it : collectEntries(dir).entrySet() )
            digestEntry(md, it.key, it.value)
        return HexFormat.of().formatHex(md.digest())
    }

    /**
     * Collect the entries that make up a plugin's executable content, keyed by relative path (with
     * forward slashes) so iteration is sorted and platform-independent. Regular files and symbolic
     * links are included: a symlink is not dereferenced here, but pf4j loads what it points to (eg.
     * a jar dropped under lib/ as a link), so it must contribute to the digest.
     */
    private static SortedMap<String,Path> collectEntries(Path dir) {
        final files = new TreeMap<String,Path>()
        Files.walkFileTree(dir, new SimpleFileVisitor<Path>() {
            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
                if( attrs.isRegularFile() || attrs.isSymbolicLink() )
                    files.put(dir.relativize(file).toString().replace('\\', '/'), file)
                return FileVisitResult.CONTINUE
            }
        })
        return files
    }

    /**
     * Feed one entry into the digest: its relative path, a type tag domain-separating a link ('L')
     * from a regular file ('F'), then the link target (not followed, so a dangling or directory
     * link neither throws nor recurses) or the file bytes.
     */
    private static void digestEntry(MessageDigest md, String rel, Path file) {
        md.update(rel.getBytes('UTF-8'))
        md.update((byte) 0)
        if( Files.isSymbolicLink(file) ) {
            md.update((byte) 'L'.charAt(0))
            md.update(Files.readSymbolicLink(file).toString().replace('\\', '/').getBytes('UTF-8'))
        }
        else {
            md.update((byte) 'F'.charAt(0))
            final buffer = new byte[8192]
            try (InputStream is = Files.newInputStream(file)) {
                int read
                while( (read = is.read(buffer)) != -1 )
                    md.update(buffer, 0, read)
            }
        }
        md.update((byte) 0)
    }
}
