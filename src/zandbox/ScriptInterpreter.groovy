/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

import groovy.transform.ToString
import groovyx.gpars.dataflow.DataflowWriteChannel

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptInterpreter {

    def List<MethodWrap> methods = []

    Binding context

    TaskProcessor processor


    @ToString(includeNames = true, includePackage = true)
    static class PropertyWrap { String name; def value;

        Object invokeMethod(String name, Object args) { new PropertyWrap(name:name, value: args) }

        Object getProperty(String propertyName) { new PropertyWrap(name: name)}

    }

    @ToString(includeNames = true, includePackage = true)
    static class MethodWrap { String name; def args }


    class PropertiesMapper {

        Object invokeMethod(String name, Object args) {
            def result = new MethodWrap(name:name, args: args)
            ScriptInterpreter.this.methods << result
            result
        }

        Object getProperty(String name) {
            def ctx = ScriptInterpreter.this.context
            def value = ctx.hasVariable(name) ? ctx.getVariable(name) : null
            def result = new PropertyWrap(name: name, value: value)
            result
        }

    }

    ScriptInterpreter( Binding binding, TaskProcessor processor ) {
        this.context = binding
        this.processor = processor
    }

    String defineChannels(Closure script) {
        pass0( script )
        pass1( script )
    }

    /**
     * First scan, detect ALL methods invocation in the black specified.
     * It populates the {@code #methods} list
     * @param closure
     * @return
     */
    private pass0(Closure closure) {
        def copy = closure.clone() as Closure
        copy.delegate = new PropertiesMapper()
        copy.call()
        methods
    }

    /**
     * Define the behavior for the methods {@code input}, {@code output} and {@code shell}
     */
    private pass1(Closure closure) {


        methods.each { MethodWrap method ->
            if ( method.name == 'input' ) {
                defineInput(method.args)
            }
            else if ( method.name == 'output' ) {
                defineOutput(method.args)
            }
        }


    }


    private void defineInput(Object args) {
        if( !args ) { throw new IllegalArgumentException("Missing attribute in 'input' script expression")}
        def argument = (args as Object[])[0]

        def name = null
        def value = null

        if( argument instanceof PropertyWrap ) {
            name = argument.@name
            value = argument.@value
        }
        else {
            throw new IllegalArgumentException("Invalid argument in 'input' declaration")
        }

        if( name && !value ) {
            throw new MissingPropertyException(name,this.getClass())
        }

        processor.input( [ (name): value] )

    }

    private defineOutput(Object args) {
        if( !args ) { throw new IllegalArgumentException("Missing attribute in 'output' script expression")}
        def argument = (args as Object[])[0]

        def name = null
        def value = null
        if ( argument instanceof PropertyWrap ) {
            name = argument.@name
            value = argument.@value
        }
        else {
            throw new IllegalArgumentException("Invalid argument in 'input' declaration")
        }

        if ( !value ) {
            processor.output(name)
            argument.@value = processor.getOutput(name)
        }
        else if ( value instanceof DataflowWriteChannel ){
            processor.output( [ (name):value ])
        }
        else {
            throw new IllegalArgumentException("Value '${name}' is not a valid as 'output' parameter -- it must be a subclass of '${DataflowWriteChannel.simpleName}'")
        }


    }

    private handleShell(Object args) {
        if( !args ) { throw new IllegalArgumentException("Missing attribute in 'shell' script expression")}
        def value = (args as Object[])[0]

        def result
        if( value instanceof PropertyWrap ) {
            result = value.name
        }
        else {
            result = value.toString()
        }

        return "\${${result}}"
    }
}
