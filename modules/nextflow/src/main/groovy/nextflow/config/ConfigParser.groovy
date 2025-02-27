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
 */

package nextflow.config

import java.nio.file.Path

/**
 * Interface for Nextflow config parsers.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
interface ConfigParser {

    /**
     * Toggle whether config include statements should be ignored.
     *
     * @param value
     */
    ConfigParser setIgnoreIncludes(boolean value)

    /**
     * Toggle whether to strip secrets when rendering the config.
     *
     * @param value
     */
    ConfigParser setStripSecrets(boolean value)

    /**
     * Toggle whether to render the source code of closures.
     *
     * @param value
     */
    ConfigParser setRenderClosureAsString(boolean value)

    /**
     * Toggle whether to raise an error if a missing property is accessed.
     *
     * @param value
     */
    ConfigParser setStrict(boolean value)

    /**
     * Define variables which will be made available to the config script.
     *
     * @param vars
     */
    ConfigParser setBinding(Map vars)

    /**
     * Define pipeline parameters which will be made available to the config script.
     *
     * @param vars
     */
    ConfigParser setParams(Map vars)

    /**
     * Parse a config object from the given source.
     */
    ConfigObject parse(String text)
    ConfigObject parse(File file)
    ConfigObject parse(Path path)

    /**
     * Set the profiles that should be applied.
     */
    ConfigParser setProfiles(List<String> profiles)

    /**
     * Get the set of available profiles.
     */
    Set<String> getProfiles()

}
