/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package test

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import org.junit.rules.ExternalResource
/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TemporaryPath extends ExternalResource {

    private final Path parentFolder;
    private Path folder;

    public TemporaryPath() {
        this(null);
    }

    public TemporaryPath(Path parentFolder) {
        this.parentFolder = parentFolder;
    }

    @Override
    protected void before() throws Throwable {
        create();
    }

    @Override
    protected void after() {
        delete();
    }

    // testing purposes only

    /**
     * for testing purposes only. Do not use.
     */
    public void create() throws IOException {
        folder = createTemporaryFolderIn(parentFolder);
    }

    /**
     * Returns a new fresh file with the given name under the temporary folder.
     */
    public Path newFile(String fileName) throws IOException {
        Path file = getRoot().resolve(fileName)
        Files.createFile(file)
        return file;
    }

    /**
     * Returns a new fresh file with a random name under the temporary folder.
     */
    public Path newFile() throws IOException {
        Files.createTempFile(getRoot(),"junit", null)
    }

    /**
     * Returns a new fresh folder with the given name under the temporary
     * folder.
     */
    public Path newFolder(String folder) throws IOException {
        return newFolder([folder] as String[]);
    }

    /**
     * Returns a new fresh folder with the given name(s) under the temporary
     * folder.
     */
    public Path newFolder(String... folderNames) throws IOException {
        Path file = getRoot();
        for (int i = 0; i < folderNames.length; i++) {
            String folderName = folderNames[i];
            validateFolderName(folderName);
            file = file.resolve(folderName);
            if (!file.mkdir() && isLastElementInArray(i, folderNames)) {
                throw new IOException("a folder with the name '$folderName' already exists");
            }
        }
        return file;
    }

    /**
     * Validates if multiple path components were used while creating a folder.
     *
     * @param folderName
     *            Name of the folder being created
     */
    private void validateFolderName(String folderName) throws IOException {
        def tempFile = Paths.get(folderName)
        if (tempFile.getParent() != null) {
            String errorMsg = "Folder name cannot consist of multiple path components separated by a file separator."
            + " Please use newFolder('MyParentFolder','MyFolder') to create hierarchies of folders";
            throw new IOException(errorMsg);
        }
    }

    private boolean isLastElementInArray(int index, String[] array) {
        return index == array.length - 1;
    }

    /**
     * Returns a new fresh folder with a random name under the temporary folder.
     */
    public Path newFolder() throws IOException {
        return createTemporaryFolderIn(getRoot());
    }

    private Path createTemporaryFolderIn(Path parentFolder) throws IOException {
        (parentFolder
                ? Files.createTempDirectory(parentFolder, "junit")
                : Files.createTempDirectory("junit") )
    }

    /**
     * @return the location of this temporary folder.
     */
    public Path getRoot() {
        if (folder == null) {
            throw new IllegalStateException("the temporary folder has not yet been created");
        }
        return folder;
    }

    /**
     * Delete all files and folders under the temporary folder. Usually not
     * called directly, since it is automatically applied by the {@link org.junit.Rule}
     */
    public void delete() {
        if (folder != null) {
            folder.deleteDir()
        }
    }



}
