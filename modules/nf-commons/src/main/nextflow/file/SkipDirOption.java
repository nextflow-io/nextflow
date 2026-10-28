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

/**
 * Copy option to skip the sub-directories with the given name when copying
 * a directory. The root directory being copied is never skipped.
 *
 * This option is only understood by Nextflow directory copy helpers, therefore
 * it must be removed with {@link #strip(CopyOption...)} before the options are
 * passed to a file system provider.
 */
public final class SkipDirOption implements CopyOption {

    private final String name;

    public SkipDirOption(String name) {
        this.name = name;
    }

    /**
     * @param source The root directory being copied
     * @param dir A directory visited while copying {@code source}
     * @return {@code true} when {@code dir} should not be copied
     */
    public boolean matches(Path source, Path dir) {
        final Path fileName = dir.getFileName();
        return fileName != null && name.equals(fileName.toString()) && !dir.equals(source);
    }

    /**
     * Remove any {@link SkipDirOption} from the given options
     */
    public static CopyOption[] strip(CopyOption... options) {
        if( options == null )
            return null;
        final CopyOption[] result = Arrays.stream(options).filter(it -> !(it instanceof SkipDirOption)).toArray(CopyOption[]::new);
        return result.length == options.length ? options : result;
    }

    @Override
    public String toString() {
        return "SkipDirOption[" + name + "]";
    }
}
