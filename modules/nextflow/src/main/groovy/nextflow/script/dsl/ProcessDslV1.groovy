/*
 * Copyright 2013-2025, Seqera Labs
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

import groovy.transform.TypeChecked
import nextflow.script.BaseScript
import nextflow.script.ProcessConfig
import nextflow.script.params.CmdEvalParam
import nextflow.script.params.EachInParam
import nextflow.script.params.EnvInParam
import nextflow.script.params.EnvOutParam
import nextflow.script.params.FileInParam
import nextflow.script.params.FileOutParam
import nextflow.script.params.StdInParam
import nextflow.script.params.StdOutParam
import nextflow.script.params.TupleInParam
import nextflow.script.params.TupleOutParam
import nextflow.script.params.ValueInParam
import nextflow.script.params.ValueOutParam

/**
 * Implements the DSL for processes.
 *
 * @see nextflow.ast.NextflowDSLImpl
 * @see nextflow.script.dsl.ProcessDsl.InputDsl
 * @see nextflow.script.dsl.ProcessDsl.OutputDsl
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@TypeChecked
class ProcessDslV1 extends ProcessBuilder {

    ProcessDslV1(BaseScript ownerScript, String processName) {
        super(new ProcessConfig(ownerScript, processName))
    }

    /// INPUTS

    void _in_each( obj ) {
        new EachInParam(config).bind(obj)
    }

    void _in_env( obj ) {
        new EnvInParam(config).bind(obj)
    }

    void _in_file( obj ) {
        new FileInParam(config).bind(obj)
    }

    void _in_path( Map opts=null, obj ) {
        new FileInParam(config)
                .setPathQualifier(true)
                .setOptions(opts)
                .bind(obj)
    }

    void _in_stdin( obj = null ) {
        final result = new StdInParam(config)
        if( obj )
            result.bind(obj)
    }

    void _in_tuple( Object... obj ) {
        new TupleInParam(config).bind(obj)
    }

    void _in_val( obj ) {
        new ValueInParam(config).bind(obj)
    }

    /// OUTPUTS

    void _out_env( Object obj ) {
        new EnvOutParam(config).bind(obj)
    }

    void _out_env( Map opts, Object obj ) {
        new EnvOutParam(config)
                .setOptions(opts)
                .bind(obj)
    }

    void _out_eval(Object obj ) {
        new CmdEvalParam(config).bind(obj)
    }

    void _out_eval(Map opts, Object obj ) {
        new CmdEvalParam(config)
            .setOptions(opts)
            .bind(obj)
    }

    void _out_file( Object obj ) {
        // note: check that is a String type to avoid to force
        // the evaluation of GString object to a string
        if( obj instanceof String && obj == '-' )
            new StdOutParam(config).bind(obj)

        else
            new FileOutParam(config).bind(obj)
    }

    void _out_path( Map opts=null, Object obj ) {
        // note: check that is a String type to avoid to force
        // the evaluation of GString object to a string
        if( obj instanceof String && obj == '-' ) {
            new StdOutParam(config)
                    .setOptions(opts)
                    .bind(obj)
        }
        else {
            new FileOutParam(config)
                    .setPathQualifier(true)
                    .setOptions(opts)
                    .bind(obj)
        }
    }

    void _out_stdout( Map opts ) {
        new StdOutParam(config)
                .setOptions(opts)
                .bind('-')
    }

    void _out_stdout( obj = null ) {
        new StdOutParam(config).bind('-')
    }

    void _out_tuple( Object... obj ) {
        new TupleOutParam(config) .bind(obj)
    }

    void _out_tuple( Map opts, Object... obj ) {
        new TupleOutParam(config)
                .setOptions(opts)
                .bind(obj)
    }

    void _out_val( Object obj ) {
        new ValueOutParam(config).bind(obj)
    }

    void _out_val( Map opts, Object obj ) {
        new ValueOutParam(config)
                .setOptions(opts)
                .bind(obj)
    }

}
