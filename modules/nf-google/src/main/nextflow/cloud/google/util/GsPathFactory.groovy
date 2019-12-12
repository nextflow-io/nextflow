/*
 * Copyright 2019, Google Inc
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

package nextflow.cloud.google.util

import java.nio.file.Path

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import groovy.transform.CompileStatic
import nextflow.file.FileSystemPathFactory

/**
 * Implements FileSystemPathFactory interface for Google storage
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class GsPathFactory extends FileSystemPathFactory {

    @Override
    protected Path parseUri(String uri) {
        if( !uri.startsWith('gs://') )
            return null
        final str = uri.substring(5)
        final p = str.indexOf('/')
        return ( p==-1
                ? CloudStorageFileSystem.forBucket(str).getPath('')
                : CloudStorageFileSystem.forBucket(str.substring(0,p)).getPath(str.substring(p)) )
    }
}
