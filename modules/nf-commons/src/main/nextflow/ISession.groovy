/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow

import java.nio.file.Path

/**
 * Nextflow session interface
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface ISession {

    /**
     * The folder where tasks temporary files are stored
     */
    Path getWorkDir()

    /**
     * The folder where the main script is contained
     */
    Path getBaseDir()

    /**
     * Holds the configuration object
     */
    Map getConfig()

    /**
     * The pipeline script name (without parent path)
     */
    String getScriptName()

    /**
     * @return List of path added to application classpath at runtime
     */
    List<Path> getLibDir()

    /**
     * The unique identifier of this session
     */
    UUID getUniqueId()

    /**
     * @return The global script binding object
     */
    Binding getBinding()

    boolean isCacheable()

    boolean isResumeMode()

}