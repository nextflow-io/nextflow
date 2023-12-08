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

package nextflow.script.dsl

import nextflow.script.BaseScript
import nextflow.script.params.*

/**
 * Implements the process DSL.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessDsl extends ProcessBuilder {

    ProcessDsl(BaseScript ownerScript, String processName) {
        super(ownerScript, processName)
    }

    /// INPUTS

    InParam _in_each( Object obj ) {
        new EachInParam(config).bind(obj)
    }

    InParam _in_env( Object obj ) {
        new EnvInParam(config).bind(obj)
    }

    InParam _in_file( Object obj ) {
        new FileInParam(config).bind(obj)
    }

    InParam _in_path( Map opts=[:], Object obj ) {
        new FileInParam(config)
                .setPathQualifier(true)
                .setOptions(opts)
                .bind(obj)
    }

    InParam _in_stdin( Object obj = null ) {
        def result = new StdInParam(config)
        if( obj )
            result.bind(obj)
        result
    }

    InParam _in_tuple( Object... obj ) {
        if( obj.length < 2 )
            throw new IllegalArgumentException("Input `tuple` must define at least two elements -- Check process `$processName`")
        new TupleInParam(config).bind(obj)
    }

    InParam _in_val( Object obj ) {
        new ValueInParam(config).bind(obj)
    }

    /// OUTPUTS

    OutParam _out_env( Object obj ) {
        new EnvOutParam(config)
                .bind(obj)
    }

    OutParam _out_env( Map opts, Object obj ) {
        new EnvOutParam(config)
                .setOptions(opts)
                .bind(obj)
    }

    OutParam _out_file( Object obj ) {
        // note: check that is a String type to avoid to force
        // the evaluation of GString object to a string
        if( obj instanceof String && obj == '-' )
            new StdOutParam(config).bind(obj)
        else
            new FileOutParam(config).bind(obj)
    }

    OutParam _out_path( Map opts=null, Object obj ) {
        // note: check that is a String type to avoid to force
        // the evaluation of GString object to a string
        if( obj instanceof String && obj == '-' )
            new StdOutParam(config)
                    .setOptions(opts)
                    .bind(obj)

        else
            new FileOutParam(config)
                    .setPathQualifier(true)
                    .setOptions(opts)
                    .bind(obj)
    }

    OutParam _out_stdout( Map opts ) {
        new StdOutParam(config)
                .setOptions(opts)
                .bind('-')
    }

    OutParam _out_stdout( obj = null ) {
        def result = new StdOutParam(config).bind('-')
        if( obj )
            result.setInto(obj)
        result
    }

    OutParam _out_tuple( Object... obj ) {
        if( obj.length < 2 )
            throw new IllegalArgumentException("Output `tuple` must define at least two elements -- Check process `$processName`")
        new TupleOutParam(config)
                .bind(obj)
    }

    OutParam _out_tuple( Map opts, Object... obj ) {
        if( obj.length < 2 )
            throw new IllegalArgumentException("Output `tuple` must define at least two elements -- Check process `$processName`")
        new TupleOutParam(config)
                .setOptions(opts)
                .bind(obj)
    }

    OutParam _out_val( Object obj ) {
        new ValueOutParam(config)
                .bind(obj)
    }

    OutParam _out_val( Map opts, Object obj ) {
        new ValueOutParam(config)
                .setOptions(opts)
                .bind(obj)
    }

}
