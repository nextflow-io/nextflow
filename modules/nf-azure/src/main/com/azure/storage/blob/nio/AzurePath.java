// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

package com.azure.storage.blob.nio;

import com.azure.core.util.logging.ClientLogger;
import com.azure.storage.blob.BlobClient;
import com.azure.storage.blob.BlobContainerClient;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.FileSystem;
import java.nio.file.InvalidPathException;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.WatchEvent;
import java.nio.file.WatchKey;
import java.nio.file.WatchService;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.stream.Stream;

/**
 * An object that may be used to locate a file in a file system.
 * <p>
 * The root component, if it is present, is the first element of the path and is denoted by a {@code ':'} as the last
 * character. Hence, only one instance of {@code ':'} may appear in a path string and it may only be the last character
 * of the first element in the path. The root component is used to identify which container a path belongs to.
 * <p>
 * Constructing a syntactically valid path does not ensure a resource exists at the given path. An error will
 * not be thrown until trying to access an invalid resource, e.g. trying to access a resource that does not exist.
 * <p>
 * Path names are case sensitive.
 * <p>
 * If a resource is accessed via a relative path, it will be resolved against the default directory of the file system.
 * The default directory is as defined in the {@link AzureFileSystem} docs.
 * <p>
 * Leading and trailing separators will be stripped from each component passed to
 * {@link AzureFileSystem#getPath(String, String...)}. This has the effect of treating "foo/" as though it were simply
 * "foo".
 */
public final class AzurePath implements Path {
    private final ClientLogger logger = new ClientLogger(AzurePath.class);
    static final String ROOT_DIR_SUFFIX = ":";

    private final AzureFileSystem parentFileSystem;
    private final String pathString;

    AzurePath(AzureFileSystem parentFileSystem, String first, String... more) {
        this.parentFileSystem = parentFileSystem;

        /*
        Break all strings into their respective elements and remove empty elements. This has the effect of stripping
        any trailing, leading, or internal delimiters so there are no duplicates/empty elements when we join.
         */
        List<String> elements = new ArrayList<>(Arrays.asList(first.split(parentFileSystem.getSeparator())));
        if (more != null) {
            for (String next : more) {
                elements.addAll(Arrays.asList(next.split(parentFileSystem.getSeparator())));
            }
        }
        elements.removeIf(String::isEmpty);

        this.pathString = String.join(this.parentFileSystem.getSeparator(), elements);

        // Validate the path string by checking usage of the reserved character ROOT_DIR_SUFFIX.
        for (int i = 0; i < elements.size(); i++) {
            String element = elements.get(i);
            /*
            If there is a root component, it must be the first element. A root component takes the format of
            "<fileStoreName>:". The ':', or ROOT_DIR_SUFFIX, if present, can only appear once, and can only be the last
            character of the first element.
             */
            if (i == 0) {
                if (element.contains(ROOT_DIR_SUFFIX) && element.indexOf(ROOT_DIR_SUFFIX) < element.length() - 1) {
                    throw LoggingUtility.logError(logger, new InvalidPathException(this.pathString, ROOT_DIR_SUFFIX
                        + " may only be used as the last character in the root component of a path"));
                }
            // No element besides the first may contain the ROOT_DIR_SUFFIX, as only the first element may be the root.
            } else if (element.contains(ROOT_DIR_SUFFIX)) {
                throw LoggingUtility.logError(logger, new InvalidPathException(this.pathString, ROOT_DIR_SUFFIX
                    + " is an invalid character except to identify the root element of this path if there is one."));
            }
        }
    }

    /**
     * Returns the file system that created this object.
     *
     * @return the file system that created this object
     */
    @Override
    public FileSystem getFileSystem() {
        return this.parentFileSystem;
    }

    /**
     * Tells whether or not this path is absolute.
     * <p>
     * An absolute path is complete in that it doesn't need to be combined with other path information in order to
     * locate a file. A path is considered absolute in this file system if it contains a root component.
     *
     * @return whether the path is absolute
     */
    @Override
    public boolean isAbsolute() {
        return this.getRoot() != null;
    }

    /**
     * Returns the root component of this path as a Path object, or null if this path does not have a root component.
     * <p>
     * The root component of this path also identifies the Azure Storage Container in which the file is stored. This
     * method will not validate that the root component corresponds to an actual file store/container in this
     * file system. It will simply return the root component of the path if one is present and syntactically valid.
     *
     * @return a path representing the root component of this path, or null
     */
    @Override
    public Path getRoot() {
        // Check if the first element of the path is formatted like a root directory.
        String[] elements = this.splitToElements();
        if (elements.length > 0 && elements[0].endsWith(ROOT_DIR_SUFFIX)) {
            return this.parentFileSystem.getPath(elements[0]);
        }
        return null;
    }

    /**
     * Returns the name of the file or directory denoted by this path as a Path object. The file name is the farthest
     * element from the root in the directory hierarchy.
     *
     * @return a path representing the name of the file or directory, or null if this path has zero elements
     */
    @Override
    public Path getFileName() {
        if (this.isRoot()) {
            return null;
        } else if (this.pathString.isEmpty()) {
            return this;
        } else {
            List<String> elements = Arrays.asList(this.splitToElements());
            return this.parentFileSystem.getPath(elements.get(elements.size() - 1));
        }
    }

    /**
     * Returns the parent path, or null if this path does not have a parent.
     * <p>
     * The parent of this path object consists of this path's root component, if any, and each element in the path
     * except for the farthest from the root in the directory hierarchy. This method does not access the file system;
     * the path or its parent may not exist. Furthermore, this method does not eliminate special names such as "." and
     * ".." that may be used in some implementations. On UNIX for example, the parent of "/a/b/c" is "/a/b", and the
     * parent of "x/y/." is "x/y". This method may be used with the normalize method, to eliminate redundant names, for
     * cases where shell-like navigation is required.
     * <p>
     * If this path has one or more elements, and no root component, then this method is equivalent to evaluating the
     * expression:
     *
     *  {@code subpath(0, getNameCount()-1);}
     *
     * @return a path representing the path's parent
     */
    @Override
    public Path getParent() {
        /*
        If this path only has one element or is empty, there is no parent. Note the root is included in the parent, so
        we don't use getNameCount here.
         */
        String[] elements = this.splitToElements();
        if (elements.length == 1 || elements.length == 0) {
            return null;
        }

        return this.parentFileSystem.getPath(
            this.pathString.substring(0, this.pathString.lastIndexOf(this.parentFileSystem.getSeparator())));
    }

    /**
     * Returns the number of name elements in the path.
     *
     * @return the number of elements in the path, or 0 if this path only represents a root component
     */
    @Override
    public int getNameCount() {
        if (this.pathString.isEmpty()) {
            return 1;
        }
        return this.splitToElements(this.withoutRoot()).length;
    }

    /**
     * Returns a name element of this path as a Path object.
     * <p>
     * The index parameter is the index of the name element to return. The element that is closest to the root in the
     * directory hierarchy has index 0. The element that is farthest from the root has index {@code count-1}.
     *
     * @param index the index of the element
     * @return the name element
     * @throws IllegalArgumentException if index is negative, index is greater than or equal to the number of elements,
     * or this path has zero name elements
     */
    @Override
    public Path getName(int index) {
        if (index < 0 || index >= this.getNameCount()) {
            throw LoggingUtility.logError(logger, new IllegalArgumentException(String.format("Index %d is out of "
                + "bounds", index)));
        }
        // If the path is empty, the only valid option is also an empty path.
        if (this.pathString.isEmpty()) {
            return this;
        }

        return this.parentFileSystem.getPath(this.splitToElements(this.withoutRoot())[index]);
    }

    /**
     * Returns a relative Path that is a subsequence of the name elements of this path.
     * <p>
     * The beginIndex and endIndex parameters specify the subsequence of name elements. The name that is closest to the
     * root in the directory hierarchy has index 0. The name that is farthest from the root has index {@code count-1}.
     * The returned Path object has the name elements that begin at beginIndex and extend to the element at index
     * {@code endIndex-1}.
     *
     * @param begin the index of the first element, inclusive
     * @param end the index of the last element, exclusive
     * @return a new Path object that is a subsequence of the name elements in this Path
     */
    @Override
    public Path subpath(int begin, int end) {
        if (begin < 0 || begin >= this.getNameCount()
            || end <= begin || end > this.getNameCount()) {
            throw LoggingUtility.logError(logger,
                new IllegalArgumentException(String.format("Values of begin: %d and end: %d are invalid", begin, end)));
        }

        String[] subnames = Stream.of(this.splitToElements(this.withoutRoot()))
            .skip(begin)
            .limit(end - begin)
            .toArray(String[]::new);

        return this.parentFileSystem.getPath(String.join(this.parentFileSystem.getSeparator(), subnames));
    }

    /**
     * Tests if this path starts with the given path.
     * <p>
     * This path starts with the given path if this path's root component starts with the root component of the given
     * path, and this path starts with the same name elements as the given path. If the given path has more name
     * elements than this path then false is returned.
     * <p>
     * If this path does not have a root component and the given path has a root component then this path does not start
     * with the given path.
     * <p>
     * If the given path is associated with a different FileSystem to this path then false is returned.
     * <p>
     * In this implementation, a root component starts with another root component if the two root components are
     * equivalent strings. In other words, if the files are stored in the same container.
     *
     * @param path the given path
     * @return true if this path starts with the given path; otherwise false
     */
    @Override
    public boolean startsWith(Path path) {
        if (!path.getFileSystem().equals(this.parentFileSystem)) {
            return false;
        }

        // An empty path never starts with another path and is never the start of another path.
        if (this.pathString.isEmpty() ^ ((AzurePath) path).pathString.isEmpty()) {
            return false;
        }

        String[] thisPathElements = this.splitToElements();
        String[] otherPathElements = ((AzurePath) path).splitToElements();
        if (otherPathElements.length > thisPathElements.length) {
            return false;
        }
        for (int i = 0; i < otherPathElements.length; i++) {
            if (!otherPathElements[i].equals(thisPathElements[i])) {
                return false;
            }
        }

        return true;
    }

    /**
     * Tests if this path starts with a Path, constructed by converting the given path string, in exactly the manner
     * specified by the startsWith(Path) method.
     *
     * @param path the given path string
     * @return true if this path starts with the given path; otherwise false
     * @throws InvalidPathException If the path string cannot be converted to a Path.
     */
    @Override
    public boolean startsWith(String path) {
        return this.startsWith(this.parentFileSystem.getPath(path));
    }

    /**
     * Tests if this path ends with the given path.
     * <p>
     * If the given path has N elements, and no root component, and this path has N or more elements, then this path
     * ends with the given path if the last N elements of each path, starting at the element farthest from the root,
     * are equal.
     * <p>
     * If the given path has a root component then this path ends with the given path if the root component of this path
     * ends with the root component of the given path, and the corresponding elements of both paths are equal. If this
     * path does not have a root component and the given path has a root component then this path does not end with the
     * given path.
     * <p>
     * If the given path is associated with a different FileSystem to this path then false is returned.
     * <p>
     * In this implementation, a root component ends with another root component if the two root components are
     * equivalent strings. In other words, if the files are stored in the same container.
     *
     * @param path the given path
     * @return true if this path ends with the given path; otherwise false
     */
    @Override
    public boolean endsWith(Path path) {
        /*
        There can only be one instance of a file system with a given id, so comparing object identity is equivalent
        to checking ids here.
         */
        if (path.getFileSystem() != this.parentFileSystem) {
            return false;
        }

        // An empty path never ends with another path and is never the end of another path.
        if (this.pathString.isEmpty() ^ ((AzurePath) path).pathString.isEmpty()) {
            return false;
        }

        String[] thisPathElements = this.splitToElements();
        String[] otherPathElements = ((AzurePath) path).splitToElements();
        if (otherPathElements.length > thisPathElements.length) {
            return false;
        }
        // If the given path has a root component, the paths must be equal.
        if (path.getRoot() != null && otherPathElements.length != thisPathElements.length) {
            return false;
        }
        for (int i = 1; i <= otherPathElements.length; i++) {
            if (!otherPathElements[otherPathElements.length - i]
                .equals(thisPathElements[thisPathElements.length - i])) {
                return false;
            }
        }
        return true;
    }

    /**
     * Tests if this path ends with a Path, constructed by converting the given path string, in exactly the manner
     * specified by the endsWith(Path) method.
     *
     * @param path the given path string
     * @return true if this path starts with the given path; otherwise false
     * @throws InvalidPathException If the path string cannot be converted to a Path.
     */
    @Override
    public boolean endsWith(String path) {
        return this.endsWith(this.parentFileSystem.getPath(path));
    }

    /**
     * Returns a path that is this path with redundant name elements eliminated.
     * <p>
     * It derives from this path, a path that does not contain redundant name elements. The "." and ".." are special
     * names used to indicate the current directory and parent directory. All occurrences of "." are considered
     * redundant. If a ".." is preceded by a non-".." name then both names are considered redundant (the process to
     * identify such names is repeated until is it no longer applicable).
     * <p>
     * This method does not access the file system; the path may not locate a file that exists. Eliminating ".." and a
     * preceding name from a path may result in the path that locates a different file than the original path
     *
     * @return the resulting path or this path if it does not contain redundant name elements; an empty path is returned
     * if this path does have a root component and all name elements are redundant
     *
     */
    @Override
    public Path normalize() {
        Deque<String> stack = new ArrayDeque<>();
        String[] pathElements = this.splitToElements();
        Path root = this.getRoot();
        String rootStr = root == null ? null : root.toString();
        for (String element : pathElements) {
            if (element.equals(".")) {
                continue;
            } else if (element.equals("..")) {
                if (rootStr != null) {
                    // Root path. We never push "..".
                    if (!stack.isEmpty() && stack.peekLast().equals(rootStr)) {
                        // Cannot go higher than root. Ignore.
                        continue;
                    } else {
                        stack.removeLast();
                    }
                } else {
                    // Relative paths can have an arbitrary number of ".." at the beginning.
                    if (stack.isEmpty()) {
                        stack.addLast(element);
                    } else if (stack.peek().equals("..")) {
                        stack.addLast(element);
                    } else {
                        stack.removeLast();
                    }
                }
            } else {
                stack.addLast(element);
            }
        }

        return this.parentFileSystem.getPath("", stack.toArray(new String[0]));
    }

    /**
     * Resolve the given path against this path.
     * <p>
     * If the other parameter is an absolute path then this method trivially returns other. If other is an empty path
     * then this method trivially returns this path. Otherwise this method considers this path to be a directory and
     * resolves the given path against this path. In the simplest case, the given path does not have a root component,
     * in which case this method joins the given path to this path and returns a resulting path that ends with the given
     * path. Where the given path has a root component then resolution is highly implementation dependent and therefore
     * unspecified.
     *
     * @param path the path to resolve against this path
     * @return the resulting path
     */
    @Override
    public Path resolve(Path path) {
        if (path.isAbsolute()) {
            return path;
        }
        if (path.getNameCount() == 0) {
            return this;
        }
        return this.parentFileSystem.getPath(this.toString(), path.toString());
    }

    /**
     * Converts a given path string to a Path and resolves it against this Path in exactly the manner specified by the
     * {@link #resolve(Path) resolve} method.
     *
     * @param path the path string to resolve against this path
     * @return the resulting path
     * @throws InvalidPathException if the path string cannot be converted to a Path.
     */
    @Override
    public Path resolve(String path) {
        return this.resolve(this.parentFileSystem.getPath(path));
    }

    /**
     * Resolves the given path against this path's parent path. This is useful where a file name needs to be replaced
     * with another file name. For example, suppose that the name separator is "/" and a path represents
     * "dir1/dir2/foo", then invoking this method with the Path "bar" will result in the Path "dir1/dir2/bar". If this
     * path does not have a parent path, or other is absolute, then this method returns other. If other is an empty path
     * then this method returns this path's parent, or where this path doesn't have a parent, the empty path.
     *
     * @param path the path to resolve against this path's parent
     * @return the resulting path
     */
    @Override
    public Path resolveSibling(Path path) {
        if (path.isAbsolute()) {
            return path;
        }

        Path parent = this.getParent();
        return parent == null ? path : parent.resolve(path);
    }

    /**
     * Converts a given path string to a Path and resolves it against this path's parent path in exactly the manner
     * specified by the resolveSibling method.
     *
     * @param path the path string to resolve against this path's parent
     * @return the resulting path
     * @throws InvalidPathException if the path string cannot be converted to a Path.
     */
    @Override
    public Path resolveSibling(String path) {
        return this.resolveSibling(this.parentFileSystem.getPath(path));
    }

    /**
     * Constructs a relative path between this path and a given path.
     * <p>
     * Relativization is the inverse of resolution. This method attempts to construct a relative path that when resolved
     * against this path, yields a path that locates the same file as the given path.
     * <p>
     * A relative path cannot be constructed if only one of the paths have a root component. If both paths have a root
     * component, it is still possible to relativize one against the other. If this path and the given path are equal
     * then an empty path is returned.
     * <p>
     * For any two normalized paths p and q, where q does not have a root component,
     *     {@code p.relativize(p.resolve(q)).equals(q)}
     *
     * @param path the path to relativize against this path
     * @return the resulting relative path, or an empty path if both paths are equal
     * @throws IllegalArgumentException if other is not a Path that can be relativized against this path
     */
    @Override
    public Path relativize(Path path) {
        if (path.getRoot() == null ^ this.getRoot() == null) {
            throw LoggingUtility.logError(logger,
                new IllegalArgumentException("Both paths must be absolute or neither can be"));
        }

        AzurePath thisNormalized = (AzurePath) this.normalize();
        Path otherNormalized = path.normalize();

        Deque<String> deque = new ArrayDeque<>(
            Arrays.asList(otherNormalized.toString().split(this.parentFileSystem.getSeparator())));

        int i = 0;
        String[] thisElements = thisNormalized.splitToElements();
        while (i < thisElements.length && !deque.isEmpty() && thisElements[i].equals(deque.peekFirst())) {
            deque.removeFirst();
            i++;
        }
        while (i < thisElements.length) {
            deque.addFirst("..");
            i++;
        }

        return this.parentFileSystem.getPath("", deque.toArray(new String[0]));
    }

    /**
     * Returns a URI to represent this path.
     * <p>
     * This method constructs an absolute URI with a scheme equal to the URI scheme that identifies the provider.
     * <p>
     * No authority component is defined for the {@code URI} returned by this method. This implementation offers the
     * same equivalence guarantee as the default provider.
     *
     * @return the URI representing this path
     * @throws SecurityException never
     */
    @Override
    public URI toUri() {
        try {
            return new URI(this.parentFileSystem.provider().getScheme(), null, "/" + this.toAbsolutePath().toString(),
                null, null);
        } catch (URISyntaxException e) {
            throw LoggingUtility.logError(logger, new IllegalStateException("Unable to create valid URI from path", e));
        }
    }

    /**
     * Returns a Path object representing the absolute path of this path.
     * <p>
     * If this path is already absolute then this method simply returns this path. Otherwise, this method resolves the
     * path against the default directory.
     *
     * @return a Path object representing the absolute path
     * @throws SecurityException never
     */
    @Override
    public Path toAbsolutePath() {
        if (this.isAbsolute()) {
            return this;
        }
        return this.parentFileSystem.getDefaultDirectory().resolve(this);
    }

    /**
     * Unsupported.
     * <p>
     * @param linkOptions options
     * @return the real path
     * @throws UnsupportedOperationException operation not suported.
     */
    @Override
    public Path toRealPath(LinkOption... linkOptions) throws IOException {
        throw new UnsupportedOperationException("Symbolic links are not supported.");
    }

    /**
     * Unsupported.
     * <p>
     * @return the file
     * @throws UnsupportedOperationException operation not suported.
     */
    @Override
    public File toFile() {
        throw new UnsupportedOperationException();
    }

    /**
     * Unsupported.
     * <p>
     * @param watchService watchService
     * @param kinds kinds
     * @param modifiers modifiers
     * @return the watch key
     * @throws UnsupportedOperationException operation not suported.
     */
    @Override
    public WatchKey register(WatchService watchService, WatchEvent.Kind<?>[] kinds, WatchEvent.Modifier... modifiers)
        throws IOException {
        throw new UnsupportedOperationException("WatchEvents are not supported.");
    }

    /**
     * Unsupported.
     * <p>
     * @param watchService watchService
     * @param kinds kinds
     * @return the watch key
     * @throws UnsupportedOperationException operation not suported.
     */
    @Override
    public WatchKey register(WatchService watchService, WatchEvent.Kind<?>... kinds) throws IOException {
        throw new UnsupportedOperationException("WatchEvents are not supported.");
    }

    /**
     * Returns an iterator over the name elements of this path.
     * <p>
     * The first element returned by the iterator represents the name element that is closest to the root in the
     * directory hierarchy, the second element is the next closest, and so on. The last element returned is the name of
     * the file or directory denoted by this path. The root component, if present, is not returned by the iterator.
     *
     * @return an iterator over the name elements of this path.
     */
    @Override
    public Iterator<Path> iterator() {
        if (this.pathString.isEmpty()) {
            return Collections.singletonList((Path) this).iterator();
        }
        return Arrays.asList(Stream.of(this.splitToElements(this.withoutRoot()))
            .map(s -> this.parentFileSystem.getPath(s))
            .toArray(Path[]::new))
            .iterator();
    }

    /**
     * Compares two abstract paths lexicographically. This method does not access the file system and neither file is
     * required to exist.
     * <p>
     * This method may not be used to compare paths that are associated with different file system providers.
     * <p>
     * This result of this method is identical to a string comparison on the underlying path strings.
     *
     * @return zero if the argument is equal to this path, a value less than zero if this path is lexicographically less
     * than the argument, or a value greater than zero if this path is lexicographically greater than the argument
     * @throws ClassCastException if the paths are associated with different providers
     */
    @Override
    public int compareTo(Path path) {
        if (!(path instanceof AzurePath)) {
            throw LoggingUtility.logError(logger, new ClassCastException("Other path is not an instance of "
                + "AzurePath."));
        }

        return this.pathString.compareTo(((AzurePath) path).pathString);
    }

    /**
     * Returns the string representation of this path.
     * <p>
     * If this path was created by converting a path string using the getPath method then the path string returned by
     * this method may differ from the original String used to create the path.
     * <p>
     * The returned path string uses the default name separator to separate names in the path.
     *
     * @return the string representation of this path
     */
    @Override
    public String toString() {
        return this.pathString;
    }

    /**
     * A path is considered equal to another path if it is associated with the same file system instance and if the
     * path strings are equivalent.
     *
     * @return true if, and only if, the given object is a Path that is identical to this Path
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        AzurePath paths = (AzurePath) o;
        return Objects.equals(parentFileSystem, paths.parentFileSystem)
            && Objects.equals(pathString, paths.pathString);
    }

    @Override
    public int hashCode() {
        return Objects.hash(parentFileSystem, pathString);
    }

    /*
    We don't store the blob client because unlike other types in this package, a Path does not actually indicate the
    existence or even validity of any remote resource. It is purely a representation of a path. Therefore, we do not
    construct the client or perform any validation until it is requested.
     */
    BlobClient toBlobClient() throws IOException {
        // Converting to an absolute path ensures there is a container to operate on even if it is the default.
        // Normalizing ensures the path is clean.
        Path root = this.normalize().toAbsolutePath().getRoot();
        if (root == null) {
            throw LoggingUtility.logError(logger,
                new IllegalStateException("Root should never be null after calling toAbsolutePath."));
        }
        String fileStoreName = this.rootToFileStore(root.toString());

        BlobContainerClient containerClient =
            ((AzureFileStore) this.parentFileSystem.getFileStore(fileStoreName)).getContainerClient();

        String blobName = this.withoutRoot();
        if (blobName.isEmpty()) {
            throw new IOException("Cannot get a blob client to a path that only contains the root or is an empty path");
        }

        return containerClient.getBlobClient(blobName);
    }

    /**
     * @return Whether this path consists of only a root component.
     */
    boolean isRoot() {
        return this.equals(this.getRoot());
    }

    private String withoutRoot() {
        Path root = this.getRoot();
        String str = this.pathString;
        if (root != null) {
            str = this.pathString.substring(root.toString().length());
        }
        if (str.startsWith(this.parentFileSystem.getSeparator())) {
            str = str.substring(1);
        }

        return str;
    }

    private String[] splitToElements() {
        return this.splitToElements(this.pathString);
    }

    private String[] splitToElements(String str) {
        String[] arr = str.split(this.parentFileSystem.getSeparator());
        /*
        This is a special case where we split after removing the root from a path that is just the root. Or otherwise
        have an empty path.
         */
        if (arr.length == 1 && arr[0].isEmpty()) {
            return new String[0];
        }
        return arr;
    }

    private String rootToFileStore(String root) {
        return root.substring(0, root.length() - 1); // Remove the ROOT_DIR_SUFFIX
    }

    static void ensureFileSystemOpen(Path p) throws IOException {
        if (!p.getFileSystem().isOpen()) {
            throw LoggingUtility.logError(((AzurePath) p).logger,
                new IOException("FileSystem for path has been closed. Path: " + p.toString()));
        }
    }
}
