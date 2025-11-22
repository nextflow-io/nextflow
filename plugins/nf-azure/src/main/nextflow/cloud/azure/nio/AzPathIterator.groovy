/*
 * Copyright 2021, Microsoft Corp
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
package nextflow.cloud.azure.nio

import java.nio.file.DirectoryStream
import java.nio.file.Path

import com.azure.storage.blob.models.BlobContainerItem
import com.azure.storage.blob.models.BlobItem
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * Iterator for Azure blob storage
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class AzPathIterator<T> implements Iterator<Path> {

    private AzFileSystem fs

    private Iterator<T> itr

    private DirectoryStream.Filter<? super Path> filter

    private AzPath next

    private AzPath origin

    AzPathIterator(AzPath origin, Iterator<T> itr, DirectoryStream.Filter<? super Path> filter) {
        this.origin = origin
        this.fs = origin.fileSystem
        this.itr = itr
        this.filter = filter
        advance()
    }

    abstract AzPath createPath(AzFileSystem fs, T item)

    private void advance() {

        AzPath result = null
        while( result == null && itr.hasNext() ) {
            def item = itr.next()
            def path = createPath(fs, item)
            if( path == origin )    // make sure to  skip the origin path
                result = null
            else if( filter )
                result = filter.accept(path) ? path : null
            else
                result = path
        }

        next = result
    }

    @Override
    boolean hasNext() {
        return next != null
    }

    @Override
    Path next() {
        def result = next
        if( result == null )
            throw new NoSuchElementException()
        advance()
        return result
    }

    @Override
    void remove() {
        throw new UnsupportedOperationException("Operation 'remove' is not supported by AzPathIterator")
    }

    /**
     * Implements a path iterator for blob object i.e. files and path in a
     * Azure cloud storage file system
     */
    static class ForBlobs extends AzPathIterator<BlobItem> {

        ForBlobs(AzPath path, Iterator<BlobItem> itr, DirectoryStream.Filter<? super Path> filter) {
            super(path, itr, filter)
        }

        @Override
        AzPath createPath(AzFileSystem fs, BlobItem item) {
            return new AzPath(fs, item)
        }
    }

    /**
     * Implements an iterator for buckets in a Google Cloud storage file system
     */
    static class ForContainers extends AzPathIterator<BlobContainerItem> {

        ForContainers(AzPath path, Iterator<BlobContainerItem> itr, DirectoryStream.Filter<? super Path> filter) {
            super(path, itr, filter)
        }

        @Override
        AzPath createPath(AzFileSystem fs, BlobContainerItem item) {
            return fs.getPath("/${item.name}").setAttributes( new AzFileAttributes(item) )
        }
    }

}
