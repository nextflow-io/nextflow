/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

package com.upplication.s3fs;

import java.io.IOException;
import java.nio.file.attribute.BasicFileAttributeView;
import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileTime;

/**
 * Implements {@link BasicFileAttributeView} for S3 file storage
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class S3FileAttributesView implements BasicFileAttributeView  {

    private S3FileAttributes target;

    S3FileAttributesView(S3FileAttributes target) {
        this.target = target;
    }

    @Override
    public String name() {
        return "basic";
    }

    @Override
    public BasicFileAttributes readAttributes() throws IOException {
        return target;
    }

    /**
     * This API is implemented is not supported but instead of throwing an exception just do nothing
     * to not break the method {@code java.nio.file.CopyMoveHelper#copyToForeignTarget(java.nio.file.Path, java.nio.file.Path, java.nio.file.CopyOption...)}
     *
     * @param lastModifiedTime
     * @param lastAccessTime
     * @param createTime
     * @throws IOException
     */
    @Override
    public void setTimes(FileTime lastModifiedTime, FileTime lastAccessTime, FileTime createTime) throws IOException {
        // not supported
    }
}
