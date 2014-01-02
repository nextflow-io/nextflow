/*
 * Copyright (c) 2012, the authors.
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

package nextflow.processor

import java.nio.file.Path

import groovy.transform.ToString
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.ast.ProcessVarRef


/**
 * 'SharedParam' Marker interface
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface SharedParam extends InParam, OutParam { }


/**
 * Model a process shared parameter that behaves both as input and output value
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString(includePackage=false, includeSuper = true)
abstract class BaseSharedParam extends BaseParam implements SharedParam {

    /** The value to which is bound this parameter in the process execution context */
    protected Object inTarget

    /** The name which is used to bind this parameter in the process execution context */
    protected String name

    /** The target channel reference */
    protected Object outTarget

    private inChannel

    private outChannel


    protected BaseSharedParam( Script script, val ) {
        super(script)

        if( val instanceof ProcessVarRef ) {
            this.name = val.name
            this.inTarget = getScriptVar(val.name)
        }
        else {
            this.inTarget = val
        }
    }

    /**
     * Lazy initializer
     *
     * @return the parameter object itself
     */
    BaseSharedParam lazyInit() {

        /*
         * initialize *input* channel
         */
        if( inTarget instanceof Closure ) {
            inChannel = sharedValToChannel( inTarget.call() )
        }
        else {
            inChannel = sharedValToChannel( inTarget )
        }

        /*
         * initialize *output* channel
         */
        if( outTarget != null ) {
            log.trace "shared output using > channel ref: $outTarget"
            outChannel = outputValToChannel( outTarget, DataflowVariable )
        }

        return this
    }

    /**
     * Implements the {@code as} keyword for the shared param declaration
     * NOTE: since {@code as} is a keyword for the groovy programming language
     * the method as to be named {@code _as}.
     *
     * A special pre-process will replace the "as" from the user script to the "_as"
     *
     * @see nextflow.ast.SourceModifierParserPlugin
     *
     * @param value
     * @return
     */
    BaseSharedParam _as( Object value ) {
        this.name = value
        return this
    }

    /**
     * Define the channel to which the parameter output has to bind
     *
     * @param value it can be a string representing a channel variable name in the script context. If
     *      the variable does not exist it creates a {@code DataflowVariable} in the script with that name.
     *      If the specified {@code value} is a {@code DataflowWriteChannel} object, use this object
     *      as the output channel
     *
     * @return The parameter object itself
     * {@see #outputValToChannel()}
     */
    BaseSharedParam to( Object value ) {
        this.outTarget = value
        return this
    }

    /**
     * The parameter name getter
     */
    String getName() { name }

    /**
     * Input channel getter
     */
    DataflowReadChannel getInChannel() {
        init()
        inChannel
    }

    /**
     * Output channel getter
     */
    DataflowWriteChannel getOutChannel() {
        init()
        outChannel
    }

}

/**
 * Represents a shared value parameter
 */

@ToString(includePackage=false, includeSuper = true)
class ValueSharedParam extends BaseSharedParam {

    ValueSharedParam ( Script script, def val ) {
        super(script, val)
    }

}


interface FileSpec {

    String getName()

    boolean getCreate()

    boolean isDirectory()

    boolean isFile()

    FileSpec create(boolean flag)

    FileSpec type( String str )

}

@ToString(includePackage=false, includeSuper = true)
class FileSharedParam extends BaseSharedParam implements FileSpec {

    private String fileName

    private boolean isVar

    FileSharedParam( Script script, Object val ) {
        super(script, val)

        // when only a string is entered, like:
        //      file 'file_something'
        //
        // it is supposed to be the file name to be used,
        // thus: clear the 'target' field and set the 'name'
        //

        if( val instanceof ProcessVarRef ) {
            isVar = true
            fileName = name
        }
        else if( val instanceof String ) {
            fileName = val
            inTarget = null
        }
        else {
            fileName = name
        }

        if( inTarget instanceof Path ) {
            fileName = inTarget.getName()
        }
    }


    FileSharedParam _as( value ) {
        if( inTarget == null && !isVar ) {
            inTarget = fileName
        }
        fileName = value
        return this
    }

    private String fType = 'file'

    private boolean fCreate = true


    String getFileName() {
        return fileName
    }

    @Override
    boolean getCreate() {
        return fCreate
    }

    @Override
    boolean isDirectory() {
        return fType == 'dir'
    }

    @Override
    boolean isFile() {
        return fType == 'file'
    }

    @Override
    FileSharedParam create(boolean flag) {
        this.fCreate = flag
        return this
    }

    @Override
    FileSharedParam type(String value) {
        assert value in ['file','dir']
        fType = value
        return this
    }

}

