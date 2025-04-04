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

package test

import java.nio.file.FileSystem
import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes

import com.google.common.jimfs.Configuration
import com.google.common.jimfs.Jimfs
import nextflow.config.control.ConfigParser
import nextflow.script.control.ScriptParser
import org.codehaus.groovy.control.messages.SyntaxErrorMessage
import org.codehaus.groovy.syntax.SyntaxException

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TestUtils {

    /**
     * Parse and analyze a script and return the list of errors.
     *
     * @param parser
     * @param contents
     */
    static List<SyntaxException> check(ScriptParser parser, String contents) {
        final source = parser.parse('main.nf', contents.stripIndent())
        parser.analyze()
        final errorCollector = source.getErrorCollector()
        if( !errorCollector.hasErrors() )
            return Collections.emptyList()
        return errorCollector.getErrors().stream()
            .filter(e -> e instanceof SyntaxErrorMessage)
            .map(e -> e.cause)
            .toList()
    }

    /**
     * Parse and analyze a config file and return the list of errors.
     *
     * @param parser
     * @param contents
     */
    static List<SyntaxException> check(ConfigParser parser, String contents) {
        final source = parser.parse('main.nf', contents.stripIndent())
        parser.analyze()
        final errorCollector = source.getErrorCollector()
        if( !errorCollector.hasErrors() )
            return Collections.emptyList()
        return errorCollector.getErrors().stream()
            .filter(e -> e instanceof SyntaxErrorMessage)
            .map(e -> e.cause)
            .toList()
    }

    private static FileSystem fs = Jimfs.newFileSystem(Configuration.unix())

    /**
     * Create an in-memory temporary directory.
     */
    static Path tempDir() {
        final tmp = fs.getPath('/tmp')
        Files.createDirectory(tmp)
        return Files.createTempDirectory(tmp, 'test')
    }

    /**
     * Create an in-memory temporary file.
     *
     * @param parent
     * @param name
     * @param contents
     */
    static Path tempFile(Path parent, String name, String contents = null) {
        final result = parent.resolve(name)
        Files.createDirectories(result.getParent())
        Files.createFile(result)
        if( contents )
            result.text = contents
        return result
    }

    static Path tempFile(String name, String contents = null) {
        return tempFile(tempDir(), name, contents)
    }

    /**
     * Delete a directory and its contents.
     *
     * @param root
     */
    static void deleteDir(Path root) {
        Files.walkFileTree(root, new SimpleFileVisitor<Path>() {
            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
                Files.delete(file)
                FileVisitResult.CONTINUE
            }
            @Override
            FileVisitResult postVisitDirectory(Path dir, IOException exc) {
                Files.delete(dir)
                FileVisitResult.CONTINUE
            }
        })
    }

}
