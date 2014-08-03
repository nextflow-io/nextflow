/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

/**
 * Application exit status
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface ExitCode {

    static final short OK = 0

    static final short MISSING_SCRIPT_FILE = 101

    static final short INVALID_COMMAND_LINE_PARAMETER = 102

    static final short SESSION_ABORTED = 103

    static final short MISSING_UNIQUE_ID = 104

    static final short MISSING_PROPERTY = 105

    static final short INVALID_CONFIG = 106

    static final short DAEMON_NOT_FOUND = 107

    static final short COMMAND_RUNTIME_ERROR = 108

    static final short UNKNOWN_ERROR = 255


}
