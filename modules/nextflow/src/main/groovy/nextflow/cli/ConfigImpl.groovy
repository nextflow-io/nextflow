/*
 * Copyright 2013-2023, Seqera Labs
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

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import nextflow.util.ConfigHelper

/**
 * CLI `config` sub-command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ConfigImpl {

    interface Options {
        String getPipeline()
        boolean getShowAllProfiles()
        String getProfile()
        boolean getPrintProperties()
        boolean getPrintFlatten()
        boolean getSort()

        ILauncherOptions getLauncherOptions()
    }

    private OutputStream stdout = System.out

    @Delegate
    private Options options

    ConfigImpl(Options options) {
        this.options = options
    }

    /* For testing purposes only */
    ConfigImpl() {}

    void run() {
        Plugins.init()
        Path base = null
        if( pipeline ) base = getBaseDir(pipeline)
        if( !base ) base = Paths.get('.')

        if( profile && showAllProfiles ) {
            throw new AbortOperationException("Option `-profile` conflicts with option `-show-profiles`")
        }

        if( printProperties && printFlatten )
            throw new AbortOperationException("Option `-flat` and `-properties` conflicts")

        final builder = new ConfigBuilder()
                .setShowClosures(true)
                .showMissingVariables(true)
                .setLauncherOptions(launcherOptions)
                .setBaseDir(base)
                .setCmdConfig(this)

        final config = builder.buildConfigObject()

        if( printProperties ) {
            printProperties0(config, stdout)
        }
        else if( printFlatten ) {
            printFlatten0(config, stdout)
        }
        else {
            printCanonical0(config, stdout)
        }

        for( String msg : builder.warnings )
            log.warn(msg)
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

        final manager = new AssetManager(path)
        manager.isLocal() ? manager.localPath.toPath() : manager.configFile?.parent

    }

}
