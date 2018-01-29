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

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.scm.AssetManager
import picocli.CommandLine

/**
 * CLI sub-command DROP
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
//@Parameters(commandDescription = "Delete the local copy of a project")
@CommandLine.Command(name = "Drop", description ="Delete the local copy of a project")
class CmdDrop extends CmdBase {

    static final public NAME = 'drop'

    //@Parameter(required=true, description = 'name of the project to drop')
    @CommandLine.Parameters(arity = "1..*", description = "name of the project to drop")    //TODO ??
    List<String> args

    //@Parameter(names='-f', description = 'Delete the repository without taking care of local changes')
    @CommandLine.Option(names=['-f'], description = 'Delete the repository without taking care of local changes')
    boolean force

    @Override
    final String getName() { NAME }

    @Override
    void run() {

        def manager = new AssetManager(args[0])
        if( !manager.localPath.exists() ) {
            throw new AbortOperationException("No match found for: ${args[0]}")
        }

        if( this.force || manager.isClean() ) {
            manager.close()
            if( !manager.localPath.deleteDir() )
                throw new AbortOperationException("Unable to delete project `${manager.project}` -- Check access permissions for path: ${manager.localPath}")
            return
        }

        throw new AbortOperationException("Local project repository contains uncommitted changes -- wont drop it")
    }
}
