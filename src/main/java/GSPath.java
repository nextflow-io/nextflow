import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.base.Predicate;
import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import javax.annotation.Nullable;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.WatchEvent;
import java.nio.file.WatchKey;
import java.nio.file.WatchService;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static com.google.common.collect.Iterables.*;
import static java.lang.String.format;

public class GSPath implements Path {

    public static final String PATH_SEPARATOR = "/";
    /**
     * bucket name
     */
    private final String bucket;
    /**
     * Parts without bucket name
     */
    private final List<String> parts;
    /**
     * actual filesystem
     */
    private GSFileSystem fileSystem;

    private GSObjectSummary objectSummary; //This is a AWS class

    public GSPath(GSFileSystem fileSystem, String path) {

        this(fileSystem, path, "");
    }

    public GSPath(GSFileSystem fileSystem, String first,
                  String ... more) {
        String bucket = null;
        List<String> parts = Lists.newArrayList(Splitter.on(PATH_SEPARATOR).split(first));

        if (first.endsWith(PATH_SEPARATOR)) {
            parts.remove(parts.size()-1);
        }

        if (first.startsWith(PATH_SEPARATOR)) { //absolute path
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
        this.parts = KeyParts.parse(keys);
        this.fileSystem = fileSystem;
    }

    private GSPath(GSFileSystem fileSystem, String bucket,
                   Iterable<String> keys) {
        this.bucket = bucket;
        this.parts = KeyParts.parse(keys);
        this.fileSystem = fileSystem;
    }

    public String getBucket() {
        return bucket;
    }

    public String getKey() {
        if (parts.isEmpty()) {
            return "";
        }

        ImmutableList.Builder<String> builder = ImmutableList
                .<String> builder().addAll(parts);

        return Joiner.on(PATH_SEPARATOR).join(builder.build());
    }

    @Override
    public GSFileSystem getFileSystem() {
        return this.fileSystem;
    }

    @Override
    public boolean isAbsolute() {
        return bucket != null;
    }

    @Override
    public Path getRoot() {
        if (isAbsolute()) {
            return new GSPath(fileSystem, bucket, ImmutableList.<String> of());
        }

        return null;
    }

    @Override
    public Path getFileName() {
        if (!parts.isEmpty()) {
            return new GSPath(fileSystem, null, parts.subList(parts.size() - 1,
                    parts.size()));
        }
        else {
            // bucket dont have fileName
            return null;
        }
    }

    @Override
    public Path getParent() {
        // bucket is not present
        if (parts.isEmpty()) {
            return null;
        }

        if (parts.size() == 1 && (bucket == null || bucket.isEmpty())){
            return null;
        }

        return new GSPath(fileSystem, bucket,
                parts.subList(0, parts.size() - 1));
    }

    @Override
    public int getNameCount() { return parts.size(); }

    @Override
    public Path getName(int index) { return new GSPath(fileSystem, null, parts.subList(index, index + 1)); }

    @Override
    public Path subpath(int beginIndex, int endIndex) {
        return new GSPath(fileSystem, null, parts.subList(beginIndex, endIndex));
    }

    @Override
    public boolean startsWith(Path other) {
        if (other.getNameCount() > this.getNameCount()){
            return false;
        }

        if (!(other instanceof GSPath)){
            return false;
        }

        GSPath path = (GSPath) other;

        if (path.parts.size() == 0 && path.bucket == null &&
                (this.parts.size() != 0 || this.bucket != null)){
            return false;
        }

        if ((path.getBucket() != null && !path.getBucket().equals(this.getBucket())) ||
                (path.getBucket() == null && this.getBucket() != null)){
            return false;
        }

        for (int i = 0; i < path.parts.size(); i++){
            if (!path.parts.get(i).equals(this.parts.get(i))){
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean startsWith(String path) {
        GSPath other = new GSPath(this.fileSystem, path);
        return this.startsWith(other);
    }

    @Override
    public boolean endsWith(Path other) {
        if (other.getNameCount() > this.getNameCount()){
            return false;
        }
        //empty
        if (other.getNameCount() > this.getNameCount()){
            return false;
        }

        if (!(other instanceof GSPath)){
            return false;
        }

        GSPath path = (GSPath) other;

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
    public boolean endsWith(String other) { return this.endsWith(new GSPath(this.fileSystem, other)); }

    @Override
    public Path normalize() { return this; }

    @Override
    public Path resolve(Path other) {
        Preconditions.checkArgument(other instanceof GSPath,
                "other must be an instance of %s", GSPath.class.getName());

        GSPath gsPath = (GSPath) other;

        if (gsPath.isAbsolute()) {
            return gsPath;
        }

        if (gsPath.parts.isEmpty()) { // other is relative and empty
            return this;
        }

        return new GSPath(fileSystem, bucket, concat(parts, gsPath.parts));
    }

    @Override
    public Path resolve(String other) { return resolve(new GSPath(this.getFileSystem(), other)); }

    @Override
    public Path resolveSibling(Path other) {
        Preconditions.checkArgument(other instanceof GSPath,
                "other must be an instance of %s", GSPath.class.getName());

        GSPath gsPath = (GSPath) other;

        Path parent = getParent();

        if (parent == null || gsPath.isAbsolute()) {
            return gsPath;
        }

        if (gsPath.parts.isEmpty()) { // other is relative and empty
            return parent;
        }

        return new GSPath(fileSystem, bucket, concat(
                parts.subList(0, parts.size() - 1), gsPath.parts));
    }

    @Override
    public Path resolveSibling(String other) { return resolveSibling(new GSPath(this.getFileSystem(), other)); }

    @Override
    public Path relativize(Path other) {
        Preconditions.checkArgument(other instanceof GSPath,
                "other must be an instance of %s", GSPath.class.getName());
        GSPath gsPath = (GSPath) other;

        if (this.equals(other)) {
            return new GSPath(this.getFileSystem(), "");
        }

        Preconditions.checkArgument(isAbsolute(),
                "Path is already relative: %s", this);
        Preconditions.checkArgument(gsPath.isAbsolute(),
                "Cannot relative against a relative path: $s", gsPath);
        Preconditions.checkArgument(bucket.equals(gsPath.getBucket()),
                "Cannot relativize paths with different buckets: '%s', '%s'",
                this, other);

        Preconditions.checkArgument(parts.size() <= gsPath.parts.size(),
                "Cannot relativeize against a parent path: '%s', '%s'",
                this, other);

        int startPart = 0;
        for (int i = 0; i < this.parts.size(); i++){
            if (this.parts.get(i).equals((gsPath.parts.get(i))){
                startPart++;
            }
        }

        List<String> resultParts = new ArrayList<>();
        for (int i = startPart; i < gsPath.parts.size(); i++){
            resultParts.add(gsPath.parts.get(i));
        }

        return new GSPath(fileSystem, null, resultParts);
    }

    @Override
    public URI toUri() {
        StringBuilder builder = new StringBuilder();
        builder.append("gs://");
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
    public File toFile() { throw new UnsupportedOperationException(); }

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
            builder.add(new GSPath(fileSystem, null, ImmutableList.of(part)));
        }

        return builder.build().iterator();
    }

    @Override
    public int compareTo(Path other) { return toString().compareTo(other.toString()); }

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

        GSPath paths = (GSPath) o;

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

}
