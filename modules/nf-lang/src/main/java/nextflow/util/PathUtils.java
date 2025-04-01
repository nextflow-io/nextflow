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
package nextflow.util;

import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.List;

public class PathUtils {

    public static boolean isPathExcluded(Path path, List<String> excludePatterns) {
        if( excludePatterns == null )
            return false;
        return excludePatterns.stream()
            .anyMatch((pattern) -> {
                if( pattern.contains("*") || pattern.contains("?") || pattern.contains("[") ) {
                    var matcher = FileSystems.getDefault().getPathMatcher("glob:" + pattern);
                    return matcher.matches(path);
                }
                else {
                    var prefix = Path.of(pattern);
                    return path.getNameCount() >= prefix.getNameCount() && prefix.equals(path.subpath(0, prefix.getNameCount()));
                }
            });
    }

}
