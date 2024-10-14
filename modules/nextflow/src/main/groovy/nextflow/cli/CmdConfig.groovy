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
import nextflow.plugin.Plugins
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

    static final List<String> FORMATS = ['flat','properties','canonical','json','yaml']

    @Parameter(description = 'project name')
    List<String> args = []

    @Parameter(names=['-a','-show-profiles'], description = 'Show all configuration profiles')
    boolean showAllProfiles

    @Parameter(names=['-profile'], description = 'Choose a configuration profile')
    String profile

    @Deprecated
    @Parameter(names = '-properties', description = 'Prints config using Java properties notation (deprecated: use `-o properties` instead)')
    boolean printProperties

    @Deprecated
    @Parameter(names = '-flat', description = 'Print config using flat notation (deprecated: use `-o flat` instead)')
    boolean printFlatten

    @Parameter(names = '-sort', description = 'Sort config attributes')
    boolean sort

    @Parameter(names = '-value', description = 'Print the value of a config option, or fail if the option is not defined')
    String printValue

    @Parameter(names = ['-o','-output'], description = 'Print the config using the specified format: canonical,properties,flat,json,yaml')
    String outputFormat

    @Override
    String getName() { NAME }

    private OutputStream stdout = System.out

    @Override
    void run() {
        Plugins.init()
        Path base = null
        if( args ) base = getBaseDir(args[0])
        if( !base ) base = Paths.get('.')

        if( profile && showAllProfiles ) {
            throw new AbortOperationException("Option `-profile` conflicts with option `-show-profiles`")
        }

        if( printProperties && printFlatten )
            throw new AbortOperationException("Option `-flat` and `-properties` conflicts each other")

        if ( printValue && printFlatten )
            throw new AbortOperationException("Option `-value` and `-flat` conflicts each other")

        if ( printValue && printProperties )
            throw new AbortOperationException("Option `-value` and `-properties` conflicts each other")

        if( printValue && outputFormat )
            throw new AbortOperationException("Option `-value` and `-output` conflicts each other")

        if( printFlatten )
            outputFormat = 'flat'

        if( printProperties )
            outputFormat = 'properties'

        final builder = new ConfigBuilder()
                .setShowClosures(true)
                .setStripSecrets(true)
                .showMissingVariables(true)
                .setOptions(launcher.options)
                .setBaseDir(base)
                .setCmdConfig(this)

        final config = builder.buildConfigObject()

        if( printValue ) {
            printValue0(config, printValue, stdout)
        }
        else if( outputFormat=='properties' ) {
            printProperties0(config, stdout)
        }
        else if( outputFormat=='flat' ) {
            printFlatten0(config, stdout)
        }
        else if( outputFormat=='yaml' ) {
            printYaml0(config, stdout)
        }
        else if( outputFormat=='json') {
            printJson0(config, stdout)
        }
        else if( !outputFormat || outputFormat=='canonical' ) {
            printCanonical0(config, stdout)
        }
        else {
            def msg = "Unknown output format: $outputFormat"
            def suggest = FORMATS.closest(outputFormat)
            if( suggest )
                msg += " - did you mean '${suggest.first()}' instead?"
            throw new AbortOperationException(msg)
        }

        for( String msg : builder.warnings )
            log.warn(msg)
    }

    /**
     * Prints a {@link ConfigObject} using Java {@link Properties} in canonical format
     * ie. any nested config object is printed within curly brackets
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
     * Prints a property of a {@link ConfigObject}.
     *
     * @param config The {@link ConfigObject} representing the parsed workflow configuration
     * @param name The {@link String} representing the property name using dot notation
     * @param output The stream where output the formatted configuration notation
     */
    @PackageScope void printValue0(ConfigObject config, String name, OutputStream output) {
        final map = config.flatten()
        if( !map.containsKey(name) )
            throw new AbortOperationException("Configuration option '$name' not found")

        output << map.get(name).toString() << '\n'
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

    @PackageScope void printJson0(ConfigObject config, OutputStream output) {
        output << ConfigHelper.toJsonString(config, sort) << '\n'
    }

    @PackageScope void printYaml0(ConfigObject config, OutputStream output) {
        output << ConfigHelper.toYamlString(config, sort)
    }

    Path getBaseDir(String path) {

        def file = Paths.get(path)
        if( file.isDirectory() )
            return file

        if( file.exists() ) {
            return file.parent ?: Paths.get('/')
        }

        final manager = new AssetManager(path)
        manager.isLocal() ? manager.localPath.toPath() : manager.configFile?.parent

    }

}
