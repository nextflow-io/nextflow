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

package nextflow.module.spi;

import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Fallback implementation of RemoteModuleResolver that is used when no other
 * implementation is found via the SPI mechanism.
 *
 * <p>This implementation throws an exception with a helpful error message
 * indicating that remote module resolution is not available.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class FallbackRemoteModuleResolver implements RemoteModuleResolver {

    @Override
    public Path resolve(String moduleName, Path baseDir) {
        if (!Files.exists(baseDir.resolve(moduleName))) {
            throw new IllegalStateException("Module '" + moduleName + "' not locally found at 'modules' folder - use 'nextflow install' to download module files");
        }
        return baseDir.resolve(moduleName).resolve("main.nf");
    }

    @Override
    public int getPriority() {
        return Integer.MIN_VALUE;  // Fallback has lowest possible priority
    }
}