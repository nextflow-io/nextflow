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
import java.nio.file.attribute.PosixFileAttributeView
import java.nio.file.attribute.PosixFileAttributes
import java.nio.file.attribute.PosixFilePermission
import java.nio.file.attribute.PosixFilePermissions
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.atomic.AtomicBoolean

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
/**
 * Validate that a plugin directory can be trusted before its content is loaded and
 * executed inside the Nextflow process.
 *
 * The plugin cache ({@code NXF_PLUGINS_DIR}, {@code $NXF_HOME/plugins}) is expected to be a
 * per-user, private trust boundary. When it is shared and writable by other users, a lower-trust
 * user can pre-populate an official pinned plugin coordinate (e.g. {@code nf-amazon-3.10.0}) with
 * attacker-controlled code that Nextflow would then load from the cache without any origin check.
 *
 * A directory is considered untrusted when a principal other than the current user can write to
 * it, that is:
 * <ul>
 *   <li>it is owned by a user other than the one running Nextflow (compared by UID), unless the
 *       owner is {@code root} or listed in {@code NXF_PLUGINS_TRUSTED_OWNERS}; or</li>
 *   <li>it is world-writable (the {@code others-write} bit is set).</li>
 * </ul>
 *
 * This accepts the mainstream deployment shapes - a private user cache (including the
 * group-writable {@code rwxrwxr-x} created under the common umask 002 with a user-private group),
 * an admin-managed read-only cache owned by {@code root} or a configured service account - while
 * rejecting an attacker-owned directory or a world-writable cache. When running as {@code root}
 * (UID 0) the owner is always trusted (the admin's responsibility) but a world-writable cache is
 * still refused.
 *
 * Note: a self-owned cache that is group-writable to a <em>shared</em> primary group is treated as
 * trusted; this residual risk is accepted because user-private-group is the modern default and
 * flagging it would produce false positives on almost every run. World-writable is always refused.
 *
 * The behaviour is controlled by the {@code NXF_PLUGINS_STRICT_MODE} environment variable:
 * {@code warn} (default) logs a warning and continues, {@code strict} aborts the run, and
 * {@code off} disables the check. On filesystems that do not support POSIX attributes, or when the
 * owner/current UID cannot be determined, the check is a no-op (fails open).
 *
 * @author Claude <noreply@anthropic.com>
 */
@Slf4j
@CompileStatic
class PluginSecurity {

    static final String MODE_WARN = 'warn'
    static final String MODE_STRICT = 'strict'
    static final String MODE_OFF = 'off'

    static final Set<String> MODES = [MODE_WARN, MODE_STRICT, MODE_OFF] as Set

    private static final long ROOT_UID = 0

    // track dirs already reported so a multi-plugin run emits at most one warning per dir
    private static final Set<String> reported = ConcurrentHashMap.newKeySet()

    // emit the invalid-mode warning at most once per run
    private static final AtomicBoolean invalidModeWarned = new AtomicBoolean()

    /**
     * Resolve the strict mode from the {@code NXF_PLUGINS_STRICT_MODE} environment variable.
     *
     * The value is normalised (lower-cased and trimmed) and validated against the known modes.
     * An unset or blank value falls back to the default silently, while an unrecognised value
     * (e.g. a typo) falls back to the default and logs a warning (once) so it is not silently
     * ignored.
     *
     * @return one of {@code warn} (default), {@code strict} or {@code off}
     */
    static String getMode() {
        final raw = SysEnv.get('NXF_PLUGINS_STRICT_MODE')
        final mode = raw?.toLowerCase()?.trim()
        if( !mode )
            return MODE_WARN
        if( !MODES.contains(mode) ) {
            if( invalidModeWarned.compareAndSet(false, true) )
                log.warn "Invalid NXF_PLUGINS_STRICT_MODE value '${raw}' - expected one of ${MODES.join(', ')}; using default '${MODE_WARN}'"
            return MODE_WARN
        }
        return mode
    }

    /**
     * The set of owner identities trusted as plugin cache owners, in addition to the current user
     * and root, parsed from the {@code NXF_PLUGINS_TRUSTED_OWNERS} environment variable. Each
     * comma-separated entry is either a numeric UID or a user name.
     *
     * @return the trimmed, non-empty trusted-owner tokens (may be empty, never {@code null})
     */
    static Set<String> getTrustedOwners() {
        final list = SysEnv.get('NXF_PLUGINS_TRUSTED_OWNERS')
        if( !list )
            return Collections.<String>emptySet()
        return list.tokenize(',')*.trim().findAll { it } as Set<String>
    }

    /**
     * Determine whether a plugin directory should be considered untrusted.
     *
     * @param perms The POSIX permissions of the directory
     * @param ownerUid The UID of the directory owner
     * @param currentUid The UID of the user running Nextflow
     * @param trustedUids Additional owner UIDs trusted beyond the current user and root
     * @return {@code true} if the directory cannot be trusted
     */
    static boolean isUntrusted(Set<PosixFilePermission> perms, long ownerUid, long currentUid, Set<Long> trustedUids) {
        final worldWritable = perms.contains(PosixFilePermission.OTHERS_WRITE)
        // running as root: trust any owner (admin's responsibility) but still flag a world-writable cache
        if( currentUid == ROOT_UID )
            return worldWritable
        final ownerTrusted = ownerUid == currentUid || ownerUid == ROOT_UID || trustedUids.contains(ownerUid)
        if( !ownerTrusted )
            return true
        // owner is trusted; group-write is fine (user-private-group / umask 002), only world-write is a risk
        return worldWritable
    }

    /**
     * Verify that the given plugin directory can be trusted, using the mode resolved from the
     * environment.
     *
     * @param dir The plugin directory (the plugins root or a {@code $id-$version} directory)
     */
    static void checkTrustedDir(Path dir) {
        checkTrustedDir(dir, getMode())
    }

    /**
     * Verify that the given plugin directory can be trusted.
     *
     * @param dir The plugin directory (the plugins root or a {@code $id-$version} directory)
     * @param mode The strict mode: {@code warn}, {@code strict} or {@code off}
     */
    static void checkTrustedDir(Path dir, String mode) {
        if( mode == MODE_OFF || dir == null )
            return
        try {
            final view = Files.getFileAttributeView(dir, PosixFileAttributeView)
            // non-POSIX filesystem (e.g. Windows, some object stores) - nothing to check
            if( view == null )
                return
            final PosixFileAttributes attrs = view.readAttributes()
            final perms = attrs.permissions()

            // resolve owner and current UID; if either is unavailable, fail open to avoid false positives
            final ownerUid = ownerUid(dir)
            final currentUid = currentUid()
            if( ownerUid == null || currentUid == null )
                return

            final trustedUids = resolveTrustedUids(attrs.owner().name, ownerUid)
            if( isUntrusted(perms, ownerUid, currentUid, trustedUids) ) {
                final owner = attrs.owner().name
                final msg = "Plugin directory '${dir}' is not secure (owner=${owner}, mode=${PosixFilePermissions.toString(perms)}) " +
                        "- it may be modified by other users, which could allow untrusted plugin code to run in your Nextflow process"
                if( mode == MODE_STRICT )
                    throw new AbortOperationException("${msg} -- refusing to load plugins. Use a private directory (chmod 700) owned by you, add the owner to NXF_PLUGINS_TRUSTED_OWNERS, or set NXF_PLUGINS_STRICT_MODE=warn to override")
                // warn at most once per offending directory
                if( reported.add(dir.toString()) )
                    log.warn "${msg} -- set NXF_PLUGINS_STRICT_MODE=strict to refuse loading, add the owner to NXF_PLUGINS_TRUSTED_OWNERS, or use a private plugins directory owned by you"
            }
        }
        catch( AbortOperationException e ) {
            throw e
        }
        catch( Exception e ) {
            // never let a permission-inspection failure break a run
            log.debug "Unable to verify plugin directory permissions for '${dir}' - ${e.message}"
        }
    }

    /**
     * Resolve the trusted owner UIDs from {@code NXF_PLUGINS_TRUSTED_OWNERS}. Numeric entries are
     * matched by UID; a name entry matching the directory owner name resolves to the owner UID.
     */
    private static Set<Long> resolveTrustedUids(String ownerName, long ownerUid) {
        final tokens = getTrustedOwners()
        if( !tokens )
            return Collections.<Long>emptySet()
        final result = new HashSet<Long>()
        for( String tok : tokens ) {
            if( tok.isLong() )
                result.add(tok.toLong())
            else if( tok == ownerName )
                result.add(ownerUid)
        }
        return result
    }

    /**
     * @return the owner UID of the given path, or {@code null} if it cannot be determined
     */
    private static Long ownerUid(Path dir) {
        try {
            final uid = Files.getAttribute(dir, 'unix:uid')
            return uid != null ? ((Number) uid).longValue() : null
        }
        catch( Exception e ) {
            log.debug "Unable to read owner UID for '${dir}' - ${e.message}"
            return null
        }
    }

    /**
     * @return the UID of the user running Nextflow, or {@code null} if it cannot be determined
     */
    private static Long currentUid() {
        try {
            return new com.sun.security.auth.module.UnixSystem().getUid()
        }
        catch( Throwable e ) {
            log.debug "Unable to determine current user UID - ${e.message}"
            return null
        }
    }
}
