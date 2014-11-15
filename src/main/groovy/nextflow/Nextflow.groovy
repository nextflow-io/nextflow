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
import java.nio.file.Path

import groovy.io.FileType
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.file.FileHelper
import nextflow.util.ArrayTuple
/**
 * Defines the main methods imported by default in the script scope
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Nextflow {

    /**
     * Create a {@code DataflowVariable} binding it to the specified value
     *
     * @param value
     * @return
     */
    static <T> DataflowVariable<T> variable( T value = null ) {
        def result = new DataflowVariable<T>()
        if( value != null ) {
            result.bind(value)
        }
        result
    }

    /**
     * Create a {@code DataflowQueue} populating with the specified values
     * <p>
     * This 'queue' data structure can be viewed as a point-to-point (1 to 1, many to 1) communication channel.
     * It allows one or more producers send messages to one reader.
     *
     * @param values
     * @return
     */
    static <T> DataflowQueue<T> channel( Collection<T> values = null ) {

        def channel = new DataflowQueue<T>()
        if ( values )  {
            // bind e
            values.each { channel << it }

            // since queue is 'finite' close it by a poison pill
            // so the operator will stop on when all values in the queue are consumed
            // (otherwise it will wait forever for a new entry)
            channel << PoisonPill.instance
        }

        return channel
    }

    /**
     * Create a {@code DataflowQueue} populating with a single value
     * <p>
     * This 'queue' data structure can be viewed as a point-to-point (1 to 1, many to 1) communication channel.
     * It allows one or more producers send messages to one reader.
     *
     * @param item
     * @return
     */
    static <T> DataflowQueue<T> channel( T... items ) {
        return channel(items as List)
    }



    /**
     * File factory utility method.
     *
     * @param name
     * @return
     */
    static def fileNamePattern( def name ) {

        if( !name ) return null

        /*
         * expand special user home '~' character
         */
        def sName = name.toString()
        if( sName == '~' ) {
            sName = System.getProperty('user.home')
        }
        else if( sName.startsWith('~'+File.separatorChar) ) {
            sName = sName.replace('~', System.getProperty('user.home'))
        }

        /*
         * split the parent path from the file name
         */
        final path = FileHelper.asPath(sName).complete()
        def base = path.getParent()
        def filePattern = path.getFileName().toString()

        /*
         * punctual file, just return it
         */
        if( !filePattern.contains('*') && !filePattern.contains('?') ) {
            return path
        }

        /*
         * when the name contains a wildcard character, it returns the list of
         * all matching files (eventually empty)
         *
         * TODO use newDirectoryStream here and glob eventually
         */
        filePattern = filePattern.replace("?", ".?").replace("*", ".*")
        def result = new LinkedList()
        base.eachFileMatch(FileType.FILES, ~/$filePattern/ ) { result << it }
        return result

    }

    static file( def fileName ) {
        assert fileName

        if( fileName instanceof Path )
            return ((Path) fileName).complete()

        if( fileName instanceof File )
            return ((File) fileName).toPath().complete()

        // default case
        return fileNamePattern(fileName?.toString())

    }

    static files( def fileName ) {
        def result = file(fileName)
        return result instanceof List ? result : [result]
    }

    /**
     * Creates a {@link ArrayTuple} object with the given open array items
     *
     * @param args The items used to created the tuple object
     * @return An instance of {@link ArrayTuple} populated with the given argument(s)
     */
    static  ArrayTuple tuple( def value ) {
        if( !value )
            return new ArrayTuple()

        new ArrayTuple( value instanceof Collection ? (Collection)value : [value] )
    }

    /**
     * Creates a {@link ArrayTuple} object with the given open array items
     *
     * @param args The items used to created the tuple object
     * @return An instance of {@link ArrayTuple} populated with the given argument(s)
     */
    static ArrayTuple tuple( Object ... args ) {
        tuple( args as List )
    }

}
