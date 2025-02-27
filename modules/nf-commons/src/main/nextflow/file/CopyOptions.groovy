/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.file;

import java.nio.file.CopyOption;
import java.nio.file.LinkOption;
import java.nio.file.StandardCopyOption

import groovy.transform.Canonical
import groovy.transform.CompileStatic;

/**
 * Parses the arguments for a file copy operation.
 */
@CompileStatic
@Canonical
class CopyOptions {
    private boolean replaceExisting = false;
    private boolean copyAttributes = false;
    private boolean followLinks = true;

    boolean replaceExisting() { replaceExisting }
    boolean copyAttributes() { copyAttributes }
    boolean followLinks() { followLinks }

    private CopyOptions() { }

    static public CopyOptions parse(CopyOption... options) {
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
