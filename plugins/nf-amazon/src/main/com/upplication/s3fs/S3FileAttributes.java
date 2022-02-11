/*
 * Copyright 2020, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Javier Arnáiz @arnaix
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package com.upplication.s3fs;

import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileTime;

import static java.lang.String.format;

public class S3FileAttributes implements BasicFileAttributes {
	
	private final FileTime lastModifiedTime;
	private final long size;
	private final boolean directory;
	private final boolean regularFile;
	private final String key;

	public S3FileAttributes(String key, FileTime lastModifiedTime, long size,
			boolean isDirectory, boolean isRegularFile) {
		this.key = key;
		this.lastModifiedTime = lastModifiedTime;
		this.size = size;
		directory = isDirectory;
		regularFile = isRegularFile;
	}

	@Override
	public FileTime lastModifiedTime() {
		return lastModifiedTime;
	}

	@Override
	public FileTime lastAccessTime() {
		return lastModifiedTime;
	}

	@Override
	public FileTime creationTime() {
		return lastModifiedTime;
	}

	@Override
	public boolean isRegularFile() {
		return regularFile;
	}

	@Override
	public boolean isDirectory() {
		return directory;
	}

	@Override
	public boolean isSymbolicLink() {
		return false;
	}

	@Override
	public boolean isOther() {
		return false;
	}

	@Override
	public long size() {
		return size;
	}

	@Override
	public Object fileKey() {
		return key;
	}

	@Override
	public String toString() {
		return format(
				"[%s: lastModified=%s, size=%s, isDirectory=%s, isRegularFile=%s]",
				key, lastModifiedTime, size, directory, regularFile);
	}
}
