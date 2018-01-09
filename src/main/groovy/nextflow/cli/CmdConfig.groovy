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
import java.nio.file.Path
import java.nio.file.Paths

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.scm.AssetManager
import org.codehaus.groovy.runtime.InvokerHelper

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
                .setOptions(launcher.options)
                .setBaseDir(base.complete())
                .setCmdConfig(this)
                .build()
                .toConfigObject()

        if( printProperties ) {
            printProperties(config, stdout)
        }
        else if( printFlatten ) {
            printFlatten(config, stdout)
        }
        else {
            printCanonical(config, stdout)
        }
    }

    /**
     * Prints a {@link ConfigObject} using Java {@link Properties} in canonical format
     * ie. any nested config object is printed withing curly brackets
     *
     * @param config The {@link ConfigObject} representing the parsed workflow configuration
     * @param output The stream where output the formatted configuration notation
     */
    protected void printCanonical(ConfigObject config, OutputStream output) {
        def writer = new PrintWriter(output)
        canonicalFormat(writer,config,0)
        writer.flush()
    }

    private static final String TAB = '   '

    private void canonicalFormat(Writer writer, ConfigObject object, int level) {

        final keys = object.keySet().sort()

        // remove all empty config objects
        final itr = keys.iterator()
        while( itr.hasNext() ) {
            final key = itr.next()
            final value = object.get(key)
            if( value instanceof ConfigObject && value.size()==0 ) {
                itr.remove()
            }
        }

        for( int i=0; i<keys.size(); i++) {
            final key = keys[i]
            final value = object.get(key)
            if( value instanceof ConfigObject ) {
                // add an extra new-line to separate simple values from a config object
                if( level==0 && i>0 ) {
                    writer.write('\n')
                }

                writer.write(TAB*level)
                writer.write(key.toString())
                writer.write(' {\n')
                canonicalFormat(writer, value, level+1)
                writer.write(TAB*level)
                writer.write('}\n')

            }
            else {
                // add a new-line to separate simple values from a previous config object
                if( level==0 && i>0 && object.get(keys[i-1]) instanceof ConfigObject) {
                    writer.write('\n')
                }

                writer.write(TAB*level)
                writer.write(key.toString())
                writer.write(' = ')
                writer.write( InvokerHelper.inspect(value) )
                writer.write('\n')
            }
        }
    }

    /**
     * Prints a {@link ConfigObject} using Java {@link Properties} format
     *
     * @param config The {@link ConfigObject} representing the parsed workflow configuration
     * @param output The stream where output the formatted configuration notation
     */
    protected void printProperties(ConfigObject config, OutputStream output) {
        def writer = new PrintWriter(output)
        writer.write( propertiesFormat(new OrderedProperties(config.toProperties())) )
        writer.flush()
    }

    /**
     * Prints a {@link ConfigObject} using properties dot notation.
     * String values are enclosed in single quote characters.
     *
     * @param config The {@link ConfigObject} representing the parsed workflow configuration
     * @param output The stream where output the formatted configuration notation
    */
    protected void printFlatten(ConfigObject config, OutputStream output) {
        def writer = new PrintWriter(output, true)
        writer.write( flattenFormat(config) )
        writer.flush()
    }

    /**
     * Prints the {@link ConfigObject} configuration object using the default notation
     *
     * @param config The {@link ConfigObject} representing the parsed workflow configuration
     * @param output The stream where output the formatted configuration notation
     */
    protected void printDefault(ConfigObject config, OutputStream output) {
        def writer = new PrintWriter(output,true)
        config.writeTo( writer )
    }

    private String flattenFormat(ConfigObject config) {
        final props = new Properties()
        config.flatten(props)
        final ordered = new OrderedProperties(props)
        final result = new StringBuilder()
        for( String name : ordered.keys() ) {
            result << name << ' = ' << InvokerHelper.inspect(ordered.get(name))  << '\n'
        }
        result.toString()
    }

    private String propertiesFormat(Properties properties) {
        def buffer = new ByteArrayOutputStream()
        properties.store(buffer,null)
        buffer.flush()

        def result = new StringBuilder()
        for( String line : buffer.toString().readLines() ) {
            if(line.startsWith('#')) continue
            result << line << '\n'
        }
        result.toString()
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

    /**
     * Extends the basic {@link Properties} to provide the ordered enumeration of keys
     */
    static class OrderedProperties extends Properties {

        OrderedProperties() {}

        OrderedProperties( Properties properties ) {
            properties.each { key, value ->
                this.put(key,value)
            }
        }

        @Override
        Enumeration<Object> keys() {
            return new Vector<>(super.keySet().sort()).elements()
        }

    }

}
