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
 * Copyright (c) 2014 Javier Arnaiz @arnaix
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

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.WatchEvent;
import java.nio.file.WatchKey;
import java.nio.file.WatchService;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import javax.annotation.Nullable;

import com.amazonaws.services.s3.model.S3ObjectId;
import com.amazonaws.services.s3.model.S3ObjectSummary;
import com.amazonaws.services.s3.model.Tag;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.base.Predicate;
import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import nextflow.file.TagAwareFile;
import static com.google.common.collect.Iterables.concat;
import static com.google.common.collect.Iterables.filter;
import static com.google.common.collect.Iterables.transform;
import static java.lang.String.format;

public class S3Path implements Path, TagAwareFile {
	
	public static final String PATH_SEPARATOR = "/";
	/**
	 * bucket name
	 */
	private final String bucket;
	/**
	 * Parts without bucket name.
	 */
	private final List<String> parts;
	/**
	 * actual filesystem
	 */
	private S3FileSystem fileSystem;

	private S3ObjectSummary objectSummary;

	private Map<String,String> tags;

	private String contentType;

	/**
	 * path must be a string of the form "/{bucket}", "/{bucket}/{key}" or just
	 * "{key}".
	 * Examples:
	 * <ul>
	 *  <li>"/{bucket}//{value}" good, empty key paths are ignored </li>
	 * <li> "//{key}" error, missing bucket</li>
	 * <li> "/" error, missing bucket </li>
	 * </ul>
	 *
	 */
	public S3Path(S3FileSystem fileSystem, String path) {
	
		this(fileSystem, path, "");
	}

    /**
     * Build an S3Path from path segments. '/' are stripped from each segment.
     * @param first should be star with a '/' and the first element is the bucket
     * @param more directories and files
     */
	public S3Path(S3FileSystem fileSystem, String first,
                   String ... more) {

        String bucket = null;
        List<String> parts = Lists.newArrayList(Splitter.on(PATH_SEPARATOR).split(first));

        if (first.endsWith(PATH_SEPARATOR)) {
            parts.remove(parts.size()-1);
        }

        if (first.startsWith(PATH_SEPARATOR)) { // absolute path
            Preconditions.checkArgument(parts.size() >= 1,
                    "path must start with bucket name");
            Preconditions.checkArgument(!parts.get(1).isEmpty(),
                    "bucket name must be not empty");

            bucket = parts.get(1);

            if (!parts.isEmpty()) {
                parts = parts.subList(2, parts.size());
            }
        }

        if (bucket != null) {
            bucket = bucket.replace("/", "");
        }

        List<String> moreSplitted = Lists.newArrayList();

        for (String part : more){
            moreSplitted.addAll(Lists.newArrayList(Splitter.on(PATH_SEPARATOR).split(part)));
        }

        parts.addAll(moreSplitted);


		this.bucket = bucket;
		this.parts = KeyParts.parse(parts);
		this.fileSystem = fileSystem;
	}

    private S3Path(S3FileSystem fileSystem, String bucket,
                   Iterable<String> keys){
        this.bucket = bucket;
        this.parts = KeyParts.parse(keys);
        this.fileSystem = fileSystem;
    }

	
	public String getBucket() {
		return bucket;
	}
	/**
	 * key for amazon without final slash.
	 * <b>note:</b> the final slash need to be added to save a directory (Amazon s3 spec)
	 */
	public String getKey() {
		if (parts.isEmpty()) {
			return "";
		}

		ImmutableList.Builder<String> builder = ImmutableList
				.<String> builder().addAll(parts);

		return Joiner.on(PATH_SEPARATOR).join(builder.build());
	}

	public S3ObjectId toS3ObjectId() {
		return new S3ObjectId(bucket, getKey());
	}

	@Override
	public S3FileSystem getFileSystem() {
		return this.fileSystem;
	}

	@Override
	public boolean isAbsolute() {
		return bucket != null;
	}

	@Override
	public Path getRoot() {
		if (isAbsolute()) {
			return new S3Path(fileSystem, bucket, ImmutableList.<String> of());
		}

		return null;
	}

	@Override
	public Path getFileName() {
		if (!parts.isEmpty()) {
			return new S3Path(fileSystem, null, parts.subList(parts.size() - 1,
					parts.size()));
		}
        else {
            // bucket dont have fileName
            return null;
        }
	}

	@Override
	public Path getParent() {
		// bucket is not present in the parts
		if (parts.isEmpty()) {
			return null;
		}

		if (parts.size() == 1 && (bucket == null || bucket.isEmpty())){
			return null;
		}

		return new S3Path(fileSystem, bucket,
				parts.subList(0, parts.size() - 1));
	}

	@Override
	public int getNameCount() {
		return parts.size();
	}

	@Override
	public Path getName(int index) {
		return new S3Path(fileSystem, null, parts.subList(index, index + 1));
	}

	@Override
	public Path subpath(int beginIndex, int endIndex) {
		return new S3Path(fileSystem, null, parts.subList(beginIndex, endIndex));
	}

	@Override
	public boolean startsWith(Path other) {
		
		if (other.getNameCount() > this.getNameCount()){
			return false;
		}
		
		if (!(other instanceof S3Path)){
			return false;
		}
		
		S3Path path = (S3Path) other;

		if (path.parts.size() == 0 && path.bucket == null &&
				(this.parts.size() != 0 || this.bucket != null)){
			return false;
		}

		if ((path.getBucket() != null && !path.getBucket().equals(this.getBucket())) ||
				(path.getBucket() == null && this.getBucket() != null)){
			return false;
		}

		for (int i = 0; i < path.parts.size() ; i++){
			if (!path.parts.get(i).equals(this.parts.get(i))){
				return false;
			}
		}
		return true;		
	}

	@Override
	public boolean startsWith(String path) {
		S3Path other = new S3Path(this.fileSystem, path);
		return this.startsWith(other);
	}

	@Override
	public boolean endsWith(Path other) {
		if (other.getNameCount() > this.getNameCount()){
			return false;
		}
		// empty
		if (other.getNameCount() == 0 && 
				this.getNameCount() != 0){
			return false;
		}
		
		if (!(other instanceof S3Path)){
			return false;
		}

		S3Path path = (S3Path) other;

		if ((path.getBucket() != null && !path.getBucket().equals(this.getBucket())) ||
				(path.getBucket() != null && this.getBucket() == null)){
			return false;
		}
		
		// check subkeys
		
		int i = path.parts.size() - 1;
		int j = this.parts.size() - 1;
		for (; i >= 0 && j >= 0 ;){
			
			if (!path.parts.get(i).equals(this.parts.get(j))){
				return false;
			}
			i--;
			j--;
		}
		return true;
	}

	@Override
	public boolean endsWith(String other) {
		return this.endsWith(new S3Path(this.fileSystem, other));
	}

	@Override
	public Path normalize() {
		return this;
	}

	@Override
	public Path resolve(Path other) {
		Preconditions.checkArgument(other instanceof S3Path,
				"other must be an instance of %s", S3Path.class.getName());

		S3Path s3Path = (S3Path) other;

		if (s3Path.isAbsolute()) {
			return s3Path;
		}

		if (s3Path.parts.isEmpty()) { // other is relative and empty
			return this;
		}

		return new S3Path(fileSystem, bucket, concat(parts, s3Path.parts));
	}

	@Override
	public Path resolve(String other) {
		return resolve(new S3Path(this.getFileSystem(), other));
	}

	@Override
	public Path resolveSibling(Path other) {
		Preconditions.checkArgument(other instanceof S3Path,
				"other must be an instance of %s", S3Path.class.getName());

		S3Path s3Path = (S3Path) other;

		Path parent = getParent();

		if (parent == null || s3Path.isAbsolute()) {
			return s3Path;
		}

		if (s3Path.parts.isEmpty()) { // other is relative and empty
			return parent;
		}

		return new S3Path(fileSystem, bucket, concat(
				parts.subList(0, parts.size() - 1), s3Path.parts));
	}

	@Override
	public Path resolveSibling(String other) {
		return resolveSibling(new S3Path(this.getFileSystem(), other));
	}

	@Override
	public Path relativize(Path other) {
		Preconditions.checkArgument(other instanceof S3Path,
				"other must be an instance of %s", S3Path.class.getName());
		S3Path s3Path = (S3Path) other;

		if (this.equals(other)) {
			return new S3Path(this.getFileSystem(), "");
		}

		Preconditions.checkArgument(isAbsolute(),
				"Path is already relative: %s", this);
		Preconditions.checkArgument(s3Path.isAbsolute(),
				"Cannot relativize against a relative path: %s", s3Path);
		Preconditions.checkArgument(bucket.equals(s3Path.getBucket()),
				"Cannot relativize paths with different buckets: '%s', '%s'",
				this, other);
		
		Preconditions.checkArgument(parts.size() <= s3Path.parts.size(),
				"Cannot relativize against a parent path: '%s', '%s'",
				this, other);
		
		
		int startPart = 0;
		for (int i = 0; i <this.parts.size() ; i++){
			if (this.parts.get(i).equals(s3Path.parts.get(i))){
				startPart++;
			}
		}
		
		List<String> resultParts = new ArrayList<>();
		for (int i = startPart; i < s3Path.parts.size(); i++){
			resultParts.add(s3Path.parts.get(i));
		}

		return new S3Path(fileSystem, null, resultParts);
	}

	@Override
	public URI toUri() {
		StringBuilder builder = new StringBuilder();
		builder.append("s3://");
		if (fileSystem.getEndpoint() != null) {
			builder.append(fileSystem.getEndpoint());
		}
		builder.append("/");
		builder.append(bucket);
		builder.append(PATH_SEPARATOR);
		builder.append(Joiner.on(PATH_SEPARATOR).join(parts));
		return URI.create(builder.toString());
	}

	@Override
	public Path toAbsolutePath() {
		if (isAbsolute()) {
			return this;
		}

		throw new IllegalStateException(format(
				"Relative path cannot be made absolute: %s", this));
	}

	@Override
	public Path toRealPath(LinkOption... options) throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public File toFile() {
		throw new UnsupportedOperationException();
	}

	@Override
	public WatchKey register(WatchService watcher, WatchEvent.Kind<?>[] events,
			WatchEvent.Modifier... modifiers) throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public WatchKey register(WatchService watcher, WatchEvent.Kind<?>... events)
			throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public Iterator<Path> iterator() {
		ImmutableList.Builder<Path> builder = ImmutableList.builder();

		for (Iterator<String> iterator = parts.iterator(); iterator.hasNext();) {
			String part = iterator.next();
			builder.add(new S3Path(fileSystem, null, ImmutableList.of(part)));
		}

		return builder.build().iterator();
	}

	@Override
	public int compareTo(Path other) {
		return toString().compareTo(other.toString());
	}

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();

		if (isAbsolute()) {
			builder.append(PATH_SEPARATOR);
			builder.append(bucket);
			builder.append(PATH_SEPARATOR);
		}

		builder.append(Joiner.on(PATH_SEPARATOR).join(parts));

		return builder.toString();
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) {
			return true;
		}
		if (o == null || getClass() != o.getClass()) {
			return false;
		}

		S3Path paths = (S3Path) o;

		if (bucket != null ? !bucket.equals(paths.bucket)
				: paths.bucket != null) {
			return false;
		}
		if (!parts.equals(paths.parts)) {
			return false;
		}

		return true;
	}

	@Override
	public int hashCode() {
		int result = bucket != null ? bucket.hashCode() : 0;
		result = 31 * result + parts.hashCode();
		return result;
	}

	/**
	 * This method returns the cached {@link S3ObjectSummary} instance if this path has been created
	 * while iterating a directory structures by the {@link S3Iterator}.
	 * <br>
	 * After calling this method the cached object is reset, so any following method invocation will return {@code null}.
	 * This is necessary to discard the object meta-data and force to reload file attributes when required.
	 *
	 * @return The cached {@link S3ObjectSummary} for this path if any.
	 */
	public S3ObjectSummary fetchObjectSummary() {
		S3ObjectSummary result = objectSummary;
		objectSummary = null;
		return result;
	}

	// note: package scope to limit the access to this setter
	void setObjectSummary(S3ObjectSummary objectSummary) {
		this.objectSummary = objectSummary;
	}

	@Override
	public void setTags(Map<String,String> tags) {
		this.tags = tags;
	}

	@Override
	public void setContentType(String type) {
		this.contentType = type;
	}

	public List<Tag> getTagsList() {
		// nothing found, just return
		if( tags==null )
			return Collections.emptyList();
		// create a list of Tag out of the Map
		List<Tag> result = new ArrayList<>();
		for( Map.Entry<String,String> entry : tags.entrySet()) {
			result.add( new Tag(entry.getKey(), entry.getValue()) );
		}
		return result;
	}

	public String getContentType() {
		return contentType;
	}

	// ~ helpers methods

	private static Function<String, String> strip(final String ... strs) {
		return new Function<String, String>() {
			public String apply(String input) {
				String res = input;
				for (String str : strs) {
					res = res.replace(str, "");
				}
				return res;
			}
		};
	}
	
	private static Predicate<String> notEmpty() {
		return new Predicate<String>() {
			@Override
			public boolean apply(@Nullable String input) {
				return input != null && !input.isEmpty();
			}
		};
	}
	/*
	 * delete redundant "/" and empty parts
	 */
	private abstract static class KeyParts{
		
		private static ImmutableList<String> parse(List<String> parts) {
			return ImmutableList.copyOf(filter(transform(parts, strip("/")), notEmpty()));
		}
		
		private static ImmutableList<String> parse(Iterable<String> parts) {
			return ImmutableList.copyOf(filter(transform(parts, strip("/")), notEmpty()));
		}
	}
}
