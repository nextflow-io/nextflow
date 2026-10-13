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

package nextflow.file;

import java.nio.file.CopyOption;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 * Copy option to skip the sub-directories with the given names when copying
 * a directory. The root directory being copied is never skipped.
 *
 * This option is only understood by Nextflow directory copy helpers, therefore
 * it must be removed with {@link #strip(CopyOption...)} before the options are
 * passed to a file system provider.
 */
public final class SkipDirsOption implements CopyOption {

    private final Set<String> names;

    private SkipDirsOption(Set<String> names) {
        this.names = names;
    }

    public static SkipDirsOption of(String... names) {
        return new SkipDirsOption(Set.of(names));
    }

    /**
     * @param source The root directory being copied
     * @param dir A directory visited while copying {@code source}
     * @return {@code true} when {@code dir} should not be copied
     */
    public boolean matches(Path source, Path dir) {
        final Path name = dir.getFileName();
        return name != null && names.contains(name.toString()) && !dir.equals(source);
    }

    /**
     * Remove any {@link SkipDirsOption} from the given options
     */
    public static CopyOption[] strip(CopyOption... options) {
        if( options == null )
            return null;
        final List<CopyOption> result = Arrays.stream(options).filter(it -> !(it instanceof SkipDirsOption)).toList();
        return result.size() == options.length ? options : result.toArray(new CopyOption[0]);
    }

    @Override
    public String toString() {
        return "SkipDirsOption" + names;
    }
}
