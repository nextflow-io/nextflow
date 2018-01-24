/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.splitter
import java.nio.charset.Charset
import java.nio.charset.CharsetEncoder
import java.nio.file.Files
import java.nio.file.Path
import java.util.zip.GZIPOutputStream

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
/**
 * A collector strategy that creates a chunks as text files
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TextFileCollector implements CollectorStrategy, CacheableCollector, Closeable {

    final private Charset charset

    private int index

    private Writer writer

    private Path currentPath

    private boolean compress

    TextFileCollector(Path baseFile, Charset charset = Charset.defaultCharset(), boolean compress=false ) {
        assert baseFile

        this.baseFile = baseFile
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

    @Override
    void add(Object record) {

        if( record == null )
            return

        if( !writer ) {
            currentPath = getNextNameFor(baseFile, ++index)
            allPaths << currentPath
            writer = getOutputWriter(currentPath, charset, compress)
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
