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
package nextflow.agent

import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic

/**
 * Pure path-containment checks for the agent {@code filesystem} tool. Writes are
 * confined to the agent work dir; reads are allowed within the work dir or any
 * of the per-invocation readable dirs (module-output dirs). Resolves symlinks and
 * normalizes so that {@code ..} traversal and symlink targets that escape the
 * sandbox are rejected.
 *
 * <p>TOCTOU note: path resolution is not atomic. The sandbox guarantee holds only if the
 * filesystem beneath the checked paths is not concurrently mutated between this check and
 * the subsequent file operation.
 */
@CompileStatic
class SandboxGuard {

    static boolean isAllowed(Path candidate, Path workDir, Collection<Path> readableDirs, boolean write) {
        if( candidate == null || workDir == null )
            return false
        final real = realOf(candidate)
        final root = realOf(workDir)
        if( isInside(real, root) )
            return true
        if( write )
            return false
        for( final Path dir : (readableDirs ?: Collections.<Path>emptyList()) ) {
            if( dir == null )
                continue
            if( isInside(real, realOf(dir)) )
                return true
        }
        return false
    }

    /**
     * Real path of an existing target, or the normalized absolute path of the
     * nearest existing ancestor joined with the remaining (non-existent) tail —
     * so a not-yet-created write target is checked against its real parent (which
     * defeats symlink escape) rather than its literal lexical path.
     */
    private static Path realOf(Path p) {
        Path abs = p.toAbsolutePath().normalize()
        if( Files.exists(abs) )
            return abs.toRealPath()
        // walk up to the nearest existing ancestor, realpath it, re-append the tail
        Path existing = abs
        final List<String> tail = new ArrayList<String>()
        while( existing != null && !Files.exists(existing) ) {
            tail.add(0, existing.getFileName().toString())
            existing = existing.getParent()
        }
        if( existing == null )
            return abs  // safe fallback: unresolved normalized abs path won't match any real workDir
        Path base = existing.toRealPath()
        for( final String seg : tail )
            base = base.resolve(seg)
        return base.normalize()
    }

    private static boolean isInside(Path child, Path root) {
        return child.startsWith(root)
    }
}
