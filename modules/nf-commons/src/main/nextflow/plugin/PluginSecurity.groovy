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
 * A directory is considered trusted only when it is owned by the current user (or root) AND it is
 * not writable by group or others. This accepts the two legitimate deployment shapes - a private
 * user cache (e.g. {@code 0700}) and an admin-managed read-only shared cache (e.g. {@code root:root
 * 0755}) - while rejecting an attacker-owned directory or a group/world-writable cache.
 *
 * The behaviour is controlled by the {@code NXF_PLUGINS_STRICT_MODE} environment variable:
 * {@code warn} (default) logs a warning and continues, {@code strict} aborts the run, and
 * {@code off} disables the check. On filesystems that do not support POSIX attributes the check is
 * a no-op.
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

    /**
     * Resolve the strict mode from the {@code NXF_PLUGINS_STRICT_MODE} environment variable.
     *
     * The value is normalised (lower-cased and trimmed) and validated against the known modes.
     * An unset or blank value falls back to the default silently, while an unrecognised value
     * (e.g. a typo) falls back to the default and logs a warning so it is not silently ignored.
     *
     * @return one of {@code warn} (default), {@code strict} or {@code off}
     */
    static String getMode() {
        final raw = SysEnv.get('NXF_PLUGINS_STRICT_MODE')
        final mode = raw?.toLowerCase()?.trim()
        if( !mode )
            return MODE_WARN
        if( !MODES.contains(mode) ) {
            log.warn "Invalid NXF_PLUGINS_STRICT_MODE value '${raw}' - expected one of ${MODES.join(', ')}; using default '${MODE_WARN}'"
            return MODE_WARN
        }
        return mode
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
            final owner = attrs.owner().name
            final me = System.getProperty('user.name')

            final writableByOthers = perms.contains(PosixFilePermission.GROUP_WRITE) || perms.contains(PosixFilePermission.OTHERS_WRITE)
            final untrustedOwner = owner != me && owner != 'root'

            if( writableByOthers || untrustedOwner ) {
                final msg = "Plugin directory '${dir}' is not secure (owner=${owner}, mode=${PosixFilePermissions.toString(perms)}) " +
                        "- it may be modified by other users, which could allow untrusted plugin code to run in your Nextflow process"
                if( mode == MODE_STRICT )
                    throw new AbortOperationException("${msg} -- refusing to load plugins. Use a private directory (chmod 700) owned by you, or set NXF_PLUGINS_STRICT_MODE=warn to override")
                log.warn "${msg} -- set NXF_PLUGINS_STRICT_MODE=strict to refuse loading, or use a private plugins directory owned by you"
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
}
