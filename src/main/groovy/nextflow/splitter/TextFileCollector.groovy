/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.splitter

import java.nio.charset.Charset
import java.nio.charset.CharsetEncoder
import java.nio.file.Files
import java.nio.file.Path
import java.util.zip.GZIPOutputStream

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.transform.TupleConstructor
/**
 * A collector strategy that creates a chunks as text files
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TextFileCollector implements CollectorStrategy, CacheableCollector, HeaderCollector, Closeable {

    @ToString
    @EqualsAndHashCode
    @TupleConstructor
    @PackageScope
    static class CachePath {
        Path path
        HashCode hash
    }

    final private Charset charset

    private int index

    private Writer writer

    private Path currentPath

    private boolean compress

    private int count

    private String header

    TextFileCollector(CachePath base, Charset charset = Charset.defaultCharset(), boolean compress=false ) {
        assert base
        assert base.path

        this.baseFile = base.path
        this.hashCode = base.hash
        this.charset = charset
        this.compress = compress
    }

    @PackageScope
    Path getNextNameFor(Path file, int index) {
        def baseName = file.getBaseName()
        def suffix = file.getExtension()
        String fileName = suffix ? "${baseName}.${index}.${suffix}" : "${baseName}.${index}"
        if( compress )
            fileName += '.gz'
        return file.resolveSibling( fileName )
    }

    void setHeader(String value) {
        this.header = value
    }

    @Override
    void add(Object record) {

        if( record == null )
            return

        if( !writer ) {
            currentPath = getNextNameFor(baseFile, ++index)
            allPaths << currentPath
            writer = getOutputWriter(currentPath, charset, compress)
        }

        if( count++==0 && header ) {
            writer.write(header)
        }
        def str = record.toString()
        writer.write(str, 0, str.length())

    }

    protected Writer getOutputWriter(Path path, Charset charset, boolean compress) {
        if( compress ) {
            CharsetEncoder encoder = charset.newEncoder()
            def zip = new GZIPOutputStream(Files.newOutputStream(path))
            return new BufferedWriter(new OutputStreamWriter(zip, encoder))
        }
        else {
            return Files.newBufferedWriter(path,charset)
        }
    }

    @Override
    def nextChunk() {
        closeWriter()
        count = 0
        def result = currentPath
        currentPath = null
        return result
    }

    @Override
    boolean hasChunk() {
        return currentPath != null
    }

    private void closeWriter() {
        if( writer ) {
            writer.flush()
            writer.closeQuietly()
            writer = null
        }
    }

    @Override
    void close() throws IOException {
        closeWriter()
        markComplete()
    }


}
