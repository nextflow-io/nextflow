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

package nextflow.cli

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.util.HistoryFile
/**
 * CLI sub-command HISTORY
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Deprecated
@CompileStatic
class CmdHistory extends CmdBase {

    static final public NAME = 'history'

    @Override
    final String getName() { NAME }

    @Override
    void run() {
        log.info "Command `history` has been deprecated -- Use `log` instead"
        HistoryFile.DEFAULT.print()
    }
}
