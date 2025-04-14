/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.script.types.shim;

import java.util.function.Consumer;

import nextflow.script.dsl.Constant;
import nextflow.script.dsl.Description;

@Description("""
    A Path is a handle for hierarchichal paths such as local files and directories, HTTP/FTP URLs, and object storage paths (e.g. Amazon S3).

    [Read more](https://nextflow.io/docs/latest/reference/stdlib.html#path)
""")
@ShimType(java.nio.file.Path.class)
public interface Path {

    // file attributes

    @Description("""
        Returns `true` if the file exists.
    """)
    boolean exists();

    @Description("""
        Gets the file name without its extension, e.g. `/some/path/file.tar.gz` -> `file.tar`.
    """)
    String getBaseName();

    @Constant("extension")
    @Description("""
        Gets the file extension, e.g. `/some/path/file.txt` -> `txt`.
    """)
    String getExtension();

    @Constant("name")
    @Description("""
        Gets the file name, e.g. `/some/path/file.txt` -> `file.txt`.
    """)
    String getName();

    @Description("""
        Gets the file name without any extension, e.g. `/some/path/file.tar.gz` -> `file`.
    """)
    String getSimpleName();

    @Constant("parent")
    @Description("""
        Gets the file parent path, e.g. `/some/path/file.txt` -> `/some/path`.
    """)
    Path getParent();

    @Constant("scheme")
    @Description("""
        Gets the file URI scheme, e.g. `s3://some-bucket/foo.txt` -> `s3`.
    """)
    String getScheme();

    @Description("""
        Returns `true` if the file is a directory.
    """)
    boolean isDirectory();

    @Description("""
        Returns `true` if the file is empty or does not exist.
    """)
    boolean isEmpty();

    @Description("""
        Returns `true` if it is a regular file (i.e. not a directory).
    """)
    boolean isFile();

    @Description("""
        Returns `true` if the file is hidden.
    """)
    boolean isHidden();

    @Description("""
        Returns `true` if the file is a symbolic link.
    """)
    boolean isLink();

    @Description("""
        Returns the file last modified timestamp in Unix time (i.e. milliseconds since January 1, 1970).
    """)
    long lastModified();

    @Description("""
        Gets the file size in bytes.
    """)
    long size();

    @Description("""
        Gets the file path along with the protocol scheme.
    """)
    String toUriString();

    // reading

    @Description("""
        Iterates over the file, applying the specified closure to each byte.
    """)
    void eachByte(Consumer<Byte> action);

    @Description("""
        Iterates over the file, applying the specified closure to each line.
    """)
    void eachLine(Consumer<String> action);

    @Constant("bytes")
    @Description("""
        Returns the file content as a byte array.
    """)
    byte[] getBytes();

    @Constant("text")
    @Description("""
        Returns the file content as a string.
    """)
    String getText();

    @Description("""
        Reads the file line by line and returns the content as a list of strings.
    """)
    List<String> readLines();

    // writing

    @Description("""
        Appends text to a file without replacing existing content.
    """)
    void append(String text);

    @Description("""
        Writes a byte array to a file. Equivalent to setting the `bytes` property.
    """)
    void setBytes(byte[] bytes);

    @Description("""
        Writes text to a file. Equivalent to setting the `text` property.
    """)
    void setText(String text);

    @Description("""
        Writes a string to a file, replacing any existing content.
    """)
    void write(String text);

    // filesystem operations

    @Description("""
        Copies a source file or directory to a target file or directory.
    """)
    void copyTo(Path target);

    @Description("""
        Deletes the file or directory at the given path, returning `true` if the operation succeeds, and `false` otherwise.
    """)
    boolean delete();

    @Description("""
        Deletes a directory and all of its contents.
    """)
    boolean deleteDir();

    @Description("""
        Returns a file's permissions using the [symbolic notation](http://en.wikipedia.org/wiki/File_system_permissions#Symbolic_notation), e.g. `'rw-rw-r--'`.
    """)
    String getPermissions();

    @Description("""
        Returns the first-level elements (files and directories) of a directory as a list of strings.
    """)
    List<String> list();

    @Description("""
        Returns the first-level elements (files and directories) of a directory as a list of Paths.
    """)
    List<Path> listFiles();

    @Description("""
        Creates a directory at the given path, returning `true` if the directory is created successfully, and `false` otherwise.
    """)
    boolean mkdir();

    @Description("""
        Creates a directory at the given path, including any nonexistent parent directories.
    """)
    boolean mkdirs();

    @Description("""
        Creates a *filesystem link* to a given path.
    """)
    Path mklink(Map opts, String linkName);

    @Description("""
        Moves a source file or directory to a target file or directory. Follows the same semantics as `copyTo()`.
    """)
    void moveTo(Path target);

    @Description("""
        Rename a file or directory.
    """)
    boolean renameTo(String target);

    @Description("""
        Sets a file's permissions using the [symbolic notation](http://en.wikipedia.org/wiki/File_system_permissions#Symbolic_notation).
    """)
    boolean setPermissions(String permissions);

    @Description("""
        Sets a file's permissions using the [numeric notation](http://en.wikipedia.org/wiki/File_system_permissions#Numeric_notation), i.e. as three digits representing the **owner**, **group**, and **other** permissions.
    """)
    boolean setPermissions(int owner, int group, int other);

    @Description("""
        Iterates through first-level directories only.
    """)
    void eachDir(Consumer<Path> action);

    @Description("""
        Iterates through directories whose names match the given filter.
    """)
    void eachDirMatch(String nameFilter, Consumer<Path> action);

    @Description("""
        Iterates through directories depth-first (regular files are ignored).
    """)
    void eachDirRecurse(Consumer<Path> action);

    @Description("""
        Iterates through first-level files and directories.
    """)
    void eachFile(Consumer<Path> action);

    @Description("""
        Iterates through files and directories whose names match the given filter.
    """)
    void eachFileMatch(String nameFilter, Consumer<Path> action);

    @Description("""
        Iterates through files and directories depth-first.
    """)
    void eachFileRecurse(Consumer<Path> action);

    // splitting

    @Description("""
        Counts the number of records in a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file. See the {ref}`operator-splitfasta` operator for available options.
    """)
    long countFasta();

    @Description("""
        Counts the number of records in a [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file. See the {ref}`operator-splitfastq` operator for available options.
    """)
    long countFastq();

    @Description("""
        Counts the number of records in a JSON file. See the {ref}`operator-splitjson` operator for available options.
    """)
    long countJson();

    @Description("""
        Counts the number of lines in a text file. See the {ref}`operator-splittext` operator for available options.
    """)
    long countLines();

    @Description("""
        Splits a CSV file into a list of records. See the {ref}`operator-splitcsv` operator for available options.
    """)
    List<?> splitCsv();

    @Description("""
        Splits a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file into a list of records. See the {ref}`operator-splitfasta` operator for available options.
    """)
    List<?> splitFasta();

    @Description("""
        Splits a [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file into a list of records. See the {ref}`operator-splitfastq` operator for available options.
    """)
    List<?> splitFastq();

    @Description("""
        Splits a JSON file into a list of records. See the {ref}`operator-splitjson` operator for available options.
    """)
    List<?> splitJson();

    @Description("""
        Splits a text file into a list of lines. See the {ref}`operator-splittext` operator for available options.
    """)
    List<String> splitText();

}
