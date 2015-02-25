/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

package nextflow.script
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
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
@InheritConstructors
abstract class BaseSharedParam extends BaseInParam implements SharedParam {

    protected intoObject

    protected OutParam.Mode mode = BasicMode.standard

    private outChannel

    /**
     * Lazy initializer
     *
     * @return the parameter object itself
     */
    void lazyInit() {

        /*
         * initialize in channel
         */
        super.lazyInit()


        /*
         * initialize *output* channel
         */
        def importName = fromObject instanceof TokenVar ? fromObject.name : null
        def hasImported = importName && binding.hasVariable(importName)

        if( intoObject instanceof TokenVar ) {
            // make sure that the output channel name is not the same as the input channel
            if( hasImported && intoObject.name == importName )
                throw new IllegalArgumentException("shared parameter 'into' cannot be same as the 'from' value -- specified value: ${intoObject.name}")

            outChannel = outputValToChannel( intoObject.name, DataflowVariable )
        }
        else if( intoObject != null ) {
            throw new IllegalArgumentException("shared parameter 'into' requires a variable identifier to be specified -- you entered: ${intoObject} [${intoObject.class.simpleName}]")
        }

    }

    protected getScriptVar( String name ) {
        getScriptVar(name,false)
    }

    protected DataflowReadChannel inputValToChannel( def value ) {

        if( value instanceof DataflowExpression ) {
            return value
        }
        else if( value instanceof DataflowReadChannel ) {
            log.warn "Using queue channel on share 'from' declaration should be avoided -- take in consideration to change declaration for share: '$name' parameter"
            value = value.toList()
        }

        def result = new DataflowVariable()
        result.bind(value)
        result
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
    BaseSharedParam into( Object value ) {
        this.intoObject = value
        return this
    }

    /**
     * Output channel getter
     */
    DataflowWriteChannel getOutChannel() {
        init()
        outChannel
    }


    def BaseSharedParam mode( def mode ) {
        this.mode = BasicMode.parseValue(mode)
        return this
    }

    OutParam.Mode getMode() { mode }

}

/**
 * Represents a shared value parameter
 */

@InheritConstructors
class ValueSharedParam extends BaseSharedParam {  }


@InheritConstructors
class FileSharedParam extends BaseSharedParam  {

    String filePattern

    /**
     * Define the file name
     */
    FileSharedParam name( name ) {
        if( name instanceof String ) {
            filePattern = name
            return this
        }

        throw new IllegalArgumentException()
    }

    String getName() {

        if( bindObject instanceof Map ) {
            def entry = bindObject.entrySet().first()
            return entry?.key
        }

        return super.getName()

    }

    String getFilePattern() {

        if( filePattern )
            return filePattern

        if( bindObject instanceof String )
            return filePattern = bindObject

        if( bindObject instanceof Map ) {
            def entry = bindObject.entrySet().first()
            return filePattern = entry?.value
        }

        return filePattern = '*'
    }

}

