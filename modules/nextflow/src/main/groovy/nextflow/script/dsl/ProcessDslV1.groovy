/*
 * Copyright 2013-2026, Seqera Labs
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
import nextflow.script.ProcessConfigV1
import nextflow.script.ProcessDef
import nextflow.script.params.CmdEvalParam
import nextflow.script.params.EachInParam
import nextflow.script.params.EnvInParam
import nextflow.script.params.EnvOutParam
import nextflow.script.params.FileInParam
import nextflow.script.params.FileOutParam
import nextflow.script.params.InputsList
import nextflow.script.params.OutputsList
import nextflow.script.params.StdInParam
import nextflow.script.params.StdOutParam
import nextflow.script.params.TupleInParam
import nextflow.script.params.TupleOutParam
import nextflow.script.params.ValueInParam
import nextflow.script.params.ValueOutParam

/**
 * Implements the DSL for legacy processes.
 *
 * @see nextflow.ast.NextflowDSLImpl
 * @see nextflow.script.dsl.ProcessDsl.InputDslV1
 * @see nextflow.script.dsl.ProcessDsl.OutputDslV1
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@TypeChecked
class ProcessDslV1 extends ProcessBuilder {

    ProcessDslV1(BaseScript ownerScript, String processName) {
        super(new ProcessConfigV1(ownerScript, processName))
    }

    private ProcessConfigV1 configV1() {
        return (ProcessConfigV1) config
    }

    /// INPUTS

    void _in_each( obj ) {
        new EachInParam(configV1()).bind(obj)
    }

    void _in_env( obj ) {
        new EnvInParam(configV1()).bind(obj)
    }

    void _in_file( obj ) {
        new FileInParam(configV1()).bind(obj)
    }

    void _in_path( Map opts=null, obj ) {
        new FileInParam(configV1())
                .setPathQualifier(true)
                .setOptions(opts)
                .bind(obj)
    }

    void _in_stdin( obj = null ) {
        final result = new StdInParam(configV1())
        if( obj )
            result.bind(obj)
    }

    void _in_tuple( Object... obj ) {
        new TupleInParam(configV1()).bind(obj)
    }

    void _in_val( obj ) {
        new ValueInParam(configV1()).bind(obj)
    }

    /// OUTPUTS

    void _out_env( Object obj ) {
        new EnvOutParam(configV1()).bind(obj)
    }

    void _out_env( Map opts, Object obj ) {
        new EnvOutParam(configV1())
                .setOptions(opts)
                .bind(obj)
    }

    void _out_eval(Object obj ) {
        new CmdEvalParam(configV1()).bind(obj)
    }

    void _out_eval(Map opts, Object obj ) {
        new CmdEvalParam(configV1())
            .setOptions(opts)
            .bind(obj)
    }

    void _out_file( Object obj ) {
        // note: check that is a String type to avoid to force
        // the evaluation of GString object to a string
        if( obj instanceof String && obj == '-' )
            new StdOutParam(configV1()).bind(obj)

        else
            new FileOutParam(configV1()).bind(obj)
    }

    void _out_path( Map opts=null, Object obj ) {
        // note: check that is a String type to avoid to force
        // the evaluation of GString object to a string
        if( obj instanceof String && obj == '-' ) {
            new StdOutParam(configV1())
                    .setOptions(opts)
                    .bind(obj)
        }
        else {
            new FileOutParam(configV1())
                    .setPathQualifier(true)
                    .setOptions(opts)
                    .bind(obj)
        }
    }

    void _out_stdout( Map opts ) {
        new StdOutParam(configV1())
                .setOptions(opts)
                .bind('-')
    }

    void _out_stdout( obj = null ) {
        new StdOutParam(configV1()).bind('-')
    }

    void _out_tuple( Object... obj ) {
        new TupleOutParam(configV1()) .bind(obj)
    }

    void _out_tuple( Map opts, Object... obj ) {
        new TupleOutParam(configV1())
                .setOptions(opts)
                .bind(obj)
    }

    void _out_val( Object obj ) {
        new ValueOutParam(configV1()).bind(obj)
    }

    void _out_val( Map opts, Object obj ) {
        new ValueOutParam(configV1())
                .setOptions(opts)
                .bind(obj)
    }

}
