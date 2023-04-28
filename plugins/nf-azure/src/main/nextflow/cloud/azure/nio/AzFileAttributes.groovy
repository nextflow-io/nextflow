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

import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileTime
import java.time.OffsetDateTime
import java.util.concurrent.TimeUnit

import com.azure.storage.blob.BlobClient
import com.azure.storage.blob.BlobContainerClient
import com.azure.storage.blob.models.BlobContainerItem
import com.azure.storage.blob.models.BlobItem
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Implements {@link BasicFileAttributes} for Azure blob path
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AzFileAttributes implements BasicFileAttributes {

    private FileTime updateTime

    private FileTime creationTime

    private boolean directory

    private long size

    private String objectId

    static AzFileAttributes root() {
        new AzFileAttributes(size: 0, objectId: '/', directory: true)
    }

    AzFileAttributes() {}


    AzFileAttributes(BlobClient client) {
        final props = client.getProperties()
        objectId = "/${client.containerName}/${client.blobName}"
        creationTime = time(props.getCreationTime())
        updateTime = time(props.getLastModified())
        directory = client.blobName.endsWith('/')
        size = props.getBlobSize()
    }

    AzFileAttributes(String containerName, BlobItem item) {
        objectId = "/${containerName}/${item.name}"
        directory = item.name.endsWith('/')
        if( !directory ) {
            creationTime = time(item.properties.getCreationTime())
            updateTime = time(item.properties.getLastModified())
            size = item.properties.getContentLength()
        }
    }

    protected AzFileAttributes(BlobContainerClient client) {
        objectId = "/$client.blobContainerName"
        updateTime = time(client.getProperties().getLastModified())
        directory = true
    }

    protected AzFileAttributes(BlobContainerItem client) {
        objectId = "/$client.name"
        updateTime = time(client.getProperties().getLastModified())
        directory = true
    }

    protected AzFileAttributes(BlobContainerClient client, String blobName) {
        objectId = "/$client.blobContainerName/$blobName"
        directory = blobName.endsWith('/')
    }


    static protected FileTime time(Long millis) {
        millis ? FileTime.from(millis, TimeUnit.MILLISECONDS) : null
    }

    static protected FileTime time(OffsetDateTime offset) {
        time(offset.toInstant().toEpochMilli())
    }


    @Override
    FileTime lastModifiedTime() {
        updateTime
    }

    @Override
    FileTime lastAccessTime() {
        return null
    }

    @Override
    FileTime creationTime() {
        creationTime
    }

    @Override
    boolean isRegularFile() {
        return !directory
    }

    @Override
    boolean isDirectory() {
        return directory
    }

    @Override
    boolean isSymbolicLink() {
        return false
    }

    @Override
    boolean isOther() {
        return false
    }

    @Override
    long size() {
        return size
    }

    @Override
    Object fileKey() {
        return objectId
    }

    @Override
    boolean equals( Object obj ) {
        if( this.class != obj?.class ) return false
        def other = (AzFileAttributes)obj
        if( creationTime() != other.creationTime() ) return false
        if( lastModifiedTime() != other.lastModifiedTime() ) return false
        if( isRegularFile() != other.isRegularFile() ) return false
        if( size() != other.size() ) return false
        return true
    }

    @Override
    int hashCode() {
        Objects.hash( creationTime(), lastModifiedTime(), isRegularFile(), size() )
    }

}
