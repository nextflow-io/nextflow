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

package nextflow.util

import java.nio.file.FileSystems
import java.nio.file.Path
import java.util.regex.Pattern

import groovy.transform.CompileStatic

/**
 * Converts host file paths into the POSIX paths at which they are visible
 * to the Linux environment running the task scripts, e.g. a Docker container
 * launched from a Windows host.
 *
 * <p>Windows drive paths are mapped as {@code C:\some\dir -> /c/some/dir}, the
 * same convention used by Docker Desktop and Git Bash. Any other path is returned
 * unchanged, therefore on Linux and macOS hosts this is the identity function.
 */
@CompileStatic
class PosixPaths {

    private static final Pattern WINDOWS_DRIVE = ~/^([A-Za-z]):(?:[\\\/](.*))?$/

    private static final String EXTENDED_PREFIX = '\\\\?\\'

    /**
     * @return {@code true} when the string is an absolute Windows drive path, e.g. {@code C:\data}
     */
    static boolean isWindowsPath(String path) {
        return path != null && WINDOWS_DRIVE.matcher(stripExtendedPrefix(path)).matches()
    }

    /**
     * @return {@code true} when the string is a Windows UNC path, e.g. {@code \\server\share}
     */
    static boolean isUncPath(String path) {
        return path != null && path.startsWith('\\\\') && !path.startsWith(EXTENDED_PREFIX)
    }

    /**
     * Map a host path string to the corresponding POSIX path
     *
     * @param path The host path string
     * @return The POSIX path, or the path itself when it's not a Windows drive path
     */
    static String toPosix(String path) {
        if( path == null )
            return null
        if( isUncPath(path) )
            throw new IllegalArgumentException("Windows UNC paths are not supported by the task environment -- map the network share to a drive letter instead; offending path: $path")
        final matcher = WINDOWS_DRIVE.matcher(stripExtendedPrefix(path))
        if( !matcher.matches() )
            return path
        final drive = matcher.group(1).toLowerCase()
        final rest = matcher.group(2)?.replace('\\', '/')
        return rest ? "/${drive}/${rest}".toString() : "/${drive}".toString()
    }

    /**
     * Map a host {@link Path} to the corresponding POSIX path. Paths that do not
     * belong to the default (local) file system are returned as they are.
     *
     * @param path The host path
     * @return The POSIX path string
     */
    static String toPosix(Path path) {
        if( path == null )
            return null
        if( path.getFileSystem() != FileSystems.getDefault() )
            return path.toString()
        return toPosix(path.toString())
    }

    private static String stripExtendedPrefix(String path) {
        return path.startsWith(EXTENDED_PREFIX) ? path.substring(EXTENDED_PREFIX.length()) : path
    }

}
