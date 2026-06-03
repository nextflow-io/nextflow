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

import java.nio.charset.StandardCharsets
import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.util.MemoryUnit

/**
 * Decides how a module tool's file output is serialized back to the LLM.
 *
 * <p>By default a file output is returned as an opaque ABSOLUTE PATH STRING (the
 * "opaque-path contract"): the LLM gets a handle it can chain to the next tool but
 * never sees file CONTENTS. That is correct for bulk/binary data (FASTA contigs,
 * BAMs) the LLM only chains between tools.
 *
 * <p>Some tools, however, emit small, structured outputs the LLM must actually
 * REASON OVER (e.g. assembly statistics as a small {@code .json} with N50 / #contigs).
 * For those, the file CONTENTS are inlined so the LLM can read the numbers. The
 * decision is made per file output at serialization time:
 *
 * <pre>
 * ext = lowercased file extension (after the last '.' in the file name, '' if none)
 * if ext not in TEXT_EXTENSIONS         -> the absolute path String   (data/binary, chainable)
 * else if size(file) > maxBytes         -> [path: &lt;abs&gt;, note: "content not inlined: ..."]
 * else if looksBinary(file)             -> the absolute path String   (safety net)
 * else                                  -> the file UTF-8 content as a String   (the LLM reads it)
 * </pre>
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ToolOutputReader {

    /**
     * File extensions whose (small) content is inlined for the LLM to reason over.
     * Everything else stays an opaque path handle.
     */
    static final Set<String> TEXT_EXTENSIONS = (['json','tsv','csv','txt','tab','yaml','yml','log','md'] as Set).asImmutable()

    /** Number of leading bytes inspected by the binary-content sniff. */
    private static final int SNIFF_BYTES = 8192

    /**
     * Decide whether to inline the file's content or return it as an opaque path handle.
     *
     * @param file     the tool-output file
     * @param maxBytes the maximum size (in bytes) of a structured file whose content is inlined
     * @return the file UTF-8 content as a String (inlined), the absolute path String (opaque
     *         handle), or a {@code [path: ..., note: ...]} Map when an inline candidate is too large
     */
    static Object readOrHandle(Path file, long maxBytes) {
        final ext = extensionOf(file)
        // unknown / non text-like format -> opaque path handle (data/binary, chainable)
        if( !TEXT_EXTENSIONS.contains(ext) )
            return pathString(file)
        // an inline candidate that is too large -> a path handle annotated with the reason
        final size = Files.size(file)
        if( size > maxBytes )
            return [path: pathString(file), note: "content not inlined: ${new MemoryUnit(size).toString()} exceeds ${new MemoryUnit(maxBytes).toString()} limit".toString()]
        // safety net: a text-like extension that nonetheless carries binary content -> path handle
        if( looksBinary(file) )
            return pathString(file)
        // small structured text -> inline the content so the LLM can reason over it
        return new String(Files.readAllBytes(file), StandardCharsets.UTF_8)
    }

    /**
     * The lowercased file extension: the part after the LAST '.' in the FILE NAME only.
     * Returns the empty string when the file name has no '.' (e.g. {@code README}, or
     * {@code a.b/c} whose last segment {@code c} has no dot).
     */
    static String extensionOf(Path file) {
        final name = file.getFileName()?.toString() ?: ''
        final dot = name.lastIndexOf('.')
        if( dot < 0 )
            return ''
        return name.substring(dot + 1).toLowerCase()
    }

    /**
     * Whether the file looks binary: read up to the first {@value #SNIFF_BYTES} bytes and
     * return {@code true} if any byte is a NUL (0x00).
     */
    static boolean looksBinary(Path file) {
        InputStream in = null
        try {
            in = Files.newInputStream(file)
            final buf = new byte[SNIFF_BYTES]
            final n = in.read(buf)
            for( int i=0; i<n; i++ ) {
                if( buf[i] == (byte) 0 )
                    return true
            }
            return false
        }
        finally {
            in?.close()
        }
    }

    private static String pathString(Path file) {
        return file.toAbsolutePath().toString()
    }
}
