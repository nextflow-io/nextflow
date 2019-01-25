/*
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

package nextflow.file;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.CopyOption;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.FileVisitOption;
import java.nio.file.FileVisitResult;
import java.nio.file.FileVisitor;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.StandardCopyOption;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.EnumSet;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Helper class to handle copy/move files and directories
 */
class CopyMoveHelper {

    private static Logger log = LoggerFactory.getLogger(CopyMoveHelper.class);

    private CopyMoveHelper() { }

    /**
     * Parses the arguments for a file copy operation.
     */
    private static class CopyOptions {
        boolean replaceExisting = false;
        boolean copyAttributes = false;
        boolean followLinks = true;

        private CopyOptions() { }

        static CopyOptions parse(CopyOption... options) {
            CopyOptions result = new CopyOptions();
            for (CopyOption option: options) {
                if (option == StandardCopyOption.REPLACE_EXISTING) {
                    result.replaceExisting = true;
                    continue;
                }
                if (option == LinkOption.NOFOLLOW_LINKS) {
                    result.followLinks = false;
                    continue;
                }
                if (option == StandardCopyOption.COPY_ATTRIBUTES) {
                    result.copyAttributes = true;
                    continue;
                }
                if (option == null)
                    throw new NullPointerException();
                throw new UnsupportedOperationException("'" + option +
                        "' is not a recognized copy option");
            }
            return result;
        }
    }

    /**
     * Converts the given array of options for moving a file to options suitable
     * for copying the file when a move is implemented as copy + delete.
     */
    private static CopyOption[] convertMoveToCopyOptions(CopyOption... options)
            throws AtomicMoveNotSupportedException
    {
        int len = options.length;
        CopyOption[] newOptions = new CopyOption[len+1];
        for (int i=0; i<len; i++) {
            CopyOption option = options[i];
            if (option == StandardCopyOption.ATOMIC_MOVE) {
                throw new AtomicMoveNotSupportedException(null, null,
                        "Atomic move between providers is not supported");
            }
            newOptions[i] = option;
        }
        newOptions[len] = LinkOption.NOFOLLOW_LINKS;
        return newOptions;
    }


    /**
     * Copy a file to a local or foreign file system
     *
     * @param source The source file path
     * @param target The target file path
     * @param foreign When {@code true} create a copy on the remote file system
     * @param options Copy options
     * @throws IOException
     */
    private static void copyFile(Path source, Path target, boolean foreign, CopyOption... options)
            throws IOException
    {

        if( !foreign ) {
            source.getFileSystem().provider().copy(source, target, options);
            return;
        }

        try (InputStream in = Files.newInputStream(source)) {
            Files.copy(in, target);
        }
    }

    /**
     * Create a directory content copy
     *
     * @param source The source directory to copy
     * @param target The target path where directory will be copied
     * @param options Copy options
     * @throws IOException
     */
    static void copyDirectory( final Path source, final Path target, final CopyOption... options )
            throws IOException
    {

        final boolean foreign = source.getFileSystem().provider() != target.getFileSystem().provider();

        FileVisitor<Path> visitor = new SimpleFileVisitor<Path>() {

            public FileVisitResult preVisitDirectory(Path current, BasicFileAttributes attr)
                    throws IOException
            {
                // get the *delta* path against the source path
                Path rel = source.relativize(current);
                String delta = rel != null ? rel.toString() : null;
                Path newFolder = delta != null ? target.resolve(delta) : target;
                if(log.isTraceEnabled())
                log.trace("Copy DIR: $current -> " + newFolder);
                // this `copy` creates the new folder, but does not copy the contained files
                Files.createDirectory(newFolder);
                return FileVisitResult.CONTINUE;
            }

            @Override
            public FileVisitResult visitFile(Path current, BasicFileAttributes attr)
                    throws IOException
            {
                // get the *delta* path against the source path
                Path rel = source.relativize(current);
                String delta = rel != null ? rel.toString() : null;
                Path newFile = delta != null ? target.resolve(delta) : target;
                if( log.isTraceEnabled())
                    log.trace("Copy file: " + current + " -> "+newFile.toUri());
                copyFile(current, newFile, foreign, options);
                return FileVisitResult.CONTINUE;
            }

        };

        Files.walkFileTree(source, EnumSet.of(FileVisitOption.FOLLOW_LINKS), Integer.MAX_VALUE, visitor);

    }

    /**
     * Simple copy for use when source and target are associated with different
     * providers
     */
    static void copyToForeignTarget(Path source, Path target, CopyOption... options)
            throws IOException
    {
        CopyOptions opts = CopyOptions.parse(options);
        LinkOption[] linkOptions = (opts.followLinks) ? new LinkOption[0] : new LinkOption[] { LinkOption.NOFOLLOW_LINKS };

        // attributes of source file
        BasicFileAttributes attrs = Files.readAttributes(source, BasicFileAttributes.class, linkOptions);
        if (attrs.isSymbolicLink())
            throw new IOException("Copying of symbolic links not supported");

        // delete target if it exists and REPLACE_EXISTING is specified
        if (opts.replaceExisting) {
            FileHelper.deletePath(target);
        }
        else if (Files.exists(target))
            throw new FileAlreadyExistsException(target.toString());

        // create directory or copy file
        if (attrs.isDirectory()) {
            copyDirectory(source, target);
        }
        else {
            copyFile(source, target, true);
        }

    }

    /**
     * Simple move implements as copy+delete for use when source and target are
     * associated with different providers
     */
    static void moveToForeignTarget(Path source, Path target, CopyOption... options)
            throws IOException
    {
        copyToForeignTarget(source, target, convertMoveToCopyOptions(options));
        FileHelper.deletePath(source);
    }

}
