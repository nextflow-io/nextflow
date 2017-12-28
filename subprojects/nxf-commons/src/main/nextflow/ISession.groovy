/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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