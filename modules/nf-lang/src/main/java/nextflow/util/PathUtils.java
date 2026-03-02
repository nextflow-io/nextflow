/*
 * Copyright 2013-2026, Seqera Labs
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
package nextflow.util;

import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.Predicate;

public class PathUtils {

    public static void visitFiles(Path root, Predicate<Path> predicate, Consumer<Path> action) throws IOException {
        Files.walkFileTree(root, new SimpleFileVisitor<Path>() {
            @Override
            public FileVisitResult preVisitDirectory(Path path, BasicFileAttributes attrs) {
                return predicate.test(path)
                    ? FileVisitResult.CONTINUE
                    : FileVisitResult.SKIP_SUBTREE;
            }
            @Override
            public FileVisitResult visitFile(Path path, BasicFileAttributes attrs) {
                if( predicate.test(path) )
                    action.accept(path.normalize());
                return FileVisitResult.CONTINUE;
            }
        });
    }

    public static boolean isExcluded(Path path, List<String> excludePatterns) {
        if( excludePatterns == null )
            return false;
        var pathStr = path.toString();
        return excludePatterns.stream().anyMatch(pattern ->
            pathStr.equals(pattern) || pathStr.endsWith("/" + pattern)
        );
    }

}
