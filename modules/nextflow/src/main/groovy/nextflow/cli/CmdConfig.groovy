/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.cli
import java.nio.file.Path
import java.nio.file.Paths

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.scm.AssetManager
import nextflow.util.ConfigHelper
/**
 *  Prints the pipeline configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Print a project configuration")
class CmdConfig extends CmdBase {

    static final public NAME = 'config'

    @Parameter(description = 'project name')
    List<String> args = []

    @Parameter(names=['-a','-show-profiles'], description = 'Show all configuration profiles')
    boolean showAllProfiles

    @Parameter(names=['-profile'], description = 'Choose a configuration profile')
    String profile

    @Parameter(names = '-properties', description = 'Prints config using Java properties notation')
    boolean printProperties

    @Parameter(names = '-flat', description = 'Print config using flat notation')
    boolean printFlatten

    @Parameter(names = '-sort', description = 'Sort config attributes')
    boolean sort


    @Override
    String getName() { NAME }

    private OutputStream stdout = System.out

    @Override
    void run() {
        Path base = null
        if( args ) base = getBaseDir(args[0])
        if( !base ) base = Paths.get('.')

        if( profile && showAllProfiles ) {
            throw new AbortOperationException("Option `-profile` conflicts with option `-show-profiles`")
        }

        if( printProperties && printFlatten )
            throw new AbortOperationException("Option `-flat` and `-properties` conflicts")

        def config = new ConfigBuilder()
                .setShowClosures(true)
                .setOptions(launcher.options)
                .setBaseDir(base)
                .setCmdConfig(this)
                .buildConfigObject()

        if( printProperties ) {
            printProperties0(config, stdout)
        }
        else if( printFlatten ) {
            printFlatten0(config, stdout)
        }
        else {
            printCanonical0(config, stdout)
        }
    }

    /**
     * Prints a {@link ConfigObject} using Java {@link Properties} in canonical format
     * ie. any nested config object is printed withing curly brackets
     *
     * @param config The {@link ConfigObject} representing the parsed workflow configuration
     * @param output The stream where output the formatted configuration notation
     */
    @PackageScope void printCanonical0(ConfigObject config, OutputStream output) {
        output << ConfigHelper.toCanonicalString(config, sort)
    }

    /**
     * Prints a {@link ConfigObject} using Java {@link Properties} format
     *
     * @param config The {@link ConfigObject} representing the parsed workflow configuration
     * @param output The stream where output the formatted configuration notation
     */
    @PackageScope void printProperties0(ConfigObject config, OutputStream output) {
        output << ConfigHelper.toPropertiesString(config, sort)
    }

    /**
     * Prints a {@link ConfigObject} using properties dot notation.
     * String values are enclosed in single quote characters.
     *
     * @param config The {@link ConfigObject} representing the parsed workflow configuration
     * @param output The stream where output the formatted configuration notation
    */
    @PackageScope void printFlatten0(ConfigObject config, OutputStream output) {
        output << ConfigHelper.toFlattenString(config, sort)
    }

    /**
     * Prints the {@link ConfigObject} configuration object using the default notation
     *
     * @param config The {@link ConfigObject} representing the parsed workflow configuration
     * @param output The stream where output the formatted configuration notation
     */
    @PackageScope void printDefault0(ConfigObject config, OutputStream output) {
        def writer = new PrintWriter(output,true)
        config.writeTo( writer )
    }


    Path getBaseDir(String path) {

        def file = Paths.get(path)
        if( file.isDirectory() )
            return file

        if( file.exists() ) {
            return file.parent ?: Paths.get('/')
        }

        def manager = new AssetManager(path)
        manager.isLocal() ? manager.localPath.toPath() : null

    }

}
