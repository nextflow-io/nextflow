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

import com.beust.jcommander.Parameters

/**
 * Self-update command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Parameters(commandDescription = "Update nextflow runtime to the latest available version")
class CmdSelfUpdate extends CmdBase {
    @Override
    String getName() { 'self-update' }

    @Override
    void run() {
        // actually it's doing nothing, the update process is managed by the external launcher script
        // this class in only necessary to print the command line the usage output
    }
}
