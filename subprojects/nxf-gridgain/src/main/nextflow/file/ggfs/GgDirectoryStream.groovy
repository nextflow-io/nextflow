/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.file.ggfs

import java.nio.file.DirectoryIteratorException
import java.nio.file.DirectoryStream
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.gridgain.grid.ggfs.GridGgfsFile
import org.gridgain.grid.ggfs.GridGgfsPath

/**
 * Implements a DirectoryStream iterator for GridGain file system
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GgDirectoryStream implements DirectoryStream<Path> {

    private DirectoryStream.Filter<? super Path> filter;
    private GridGgfsPath path;
    private Iterator<GridGgfsFile> target;
    private GgFileSystem fileSystem

    /**
     * Create the directory stream iterator
     *
     * @param path The {@code GgPath} to iterate over
     * @param filter The {@code DirectoryStream.Filter} to be applied
     * @throws IOException
     */
    public GgDirectoryStream( GgPath path, DirectoryStream.Filter<? super Path> filter ) throws IOException {
        this.path = path.toGridGgfsPath();
        this.fileSystem = path.fileSystem
        this.filter = filter;
        this.target = path.getFileSystem().getGgfs().listFiles(path.toGridGgfsPath()).iterator()
    }


    @Override
    public Iterator<Path> iterator()  {

        return new Iterator<Path>() {

            Path nextValue = findNext(target);

            @Override
            public boolean hasNext() {
                return nextValue != null;
            }

            @Override
            public Path next() {
                Path result = nextValue;
                nextValue = findNext(target);
                return result;
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };

    }

    protected Path findNext(Iterator<GridGgfsFile> it)  {
        while( it.hasNext() ) {
            GridGgfsFile item = it.next();
            try {
                def path = new GgPath(fileSystem, item.path().toString())
                if( filter == null || filter.accept(path) ) {
                    return path;
                }
            }
            catch( IOException e ) {
                throw new DirectoryIteratorException(e);
            }
        }

        return null;
    }

    @Override
    public void close() throws IOException {
        // nothing to do
    }


}
