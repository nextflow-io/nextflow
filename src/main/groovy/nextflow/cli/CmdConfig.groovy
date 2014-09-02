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

package nextflow.cli
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.scm.AssetManager
import nextflow.script.ConfigBuilder
/**
 *  Prints the pipeline configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Print a pipeline configuration")
class CmdConfig extends CmdX {

    static final NAME = 'config'

    @Parameter(description = 'pipeline name')
    List<String> args = []

    @Override
    String getName() { NAME }

    @Override
    void run() {
        File base = args ? getBaseDir(args[0]) : null
        def config = new ConfigBuilder()
                .setOptions(launcher.options)
                .setBaseDir(base)
                .build()

        PrintWriter stdout = new PrintWriter(System.out,true);
        config.writeTo( stdout )
    }


    File getBaseDir(String path) {

        def file = new File(path)
        if( file.isDirectory() )
            return file

        if( file.exists() ) {
            file.canonicalFile.parentFile ?: new File('/')
        }

        def manager = new AssetManager(path)
        manager.isLocal() ? manager.localPath : null

    }
}
