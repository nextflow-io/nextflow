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

package nextflow.script.dsl

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import nextflow.processor.TaskOutputCollector
import nextflow.script.BaseScript
import nextflow.script.LazyAware
import nextflow.script.LazyList
import nextflow.script.LazyVar
import nextflow.script.ProcessDef
import nextflow.script.ProcessFileInput
import nextflow.script.ProcessFileOutput
import nextflow.script.ProcessInputs
import nextflow.script.ProcessOutput
import nextflow.script.ProcessOutputs
import nextflow.script.TokenEnvCall
import nextflow.script.TokenEvalCall
import nextflow.script.TokenFileCall
import nextflow.script.TokenPathCall
import nextflow.script.TokenStdinCall
import nextflow.script.TokenStdoutCall
import nextflow.script.TokenValCall

/**
 * Implements the process DSL.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessDsl extends ProcessBuilder {

    private ProcessInputs inputs = new ProcessInputs()

    private ProcessOutputs outputs = new ProcessOutputs()

    ProcessDsl(BaseScript ownerScript, String processName) {
        super(ownerScript, processName)
    }

    /// INPUTS

    void _in_each(LazyVar var) {
        _in_val(var)
        inputs.last().setIterator(true)
    }

    void _in_each(TokenFileCall file) {
        _in_file(file.target)
        inputs.last().setIterator(true)
    }

    void _in_each(TokenPathCall path) {
        _in_path(path.target)
        inputs.last().setIterator(true)
    }

    void _in_env(LazyVar var) {
        final param = "\$in${inputs.size()}".toString()
        inputs.addParam(param)
        inputs.addEnv(var.name, new LazyVar(param))
    }

    void _in_file(Object source) {
        final param = _in_path0(source, false, [:])
        inputs.addParam(param)
    }

    void _in_path(Map opts=[:], Object source) {
        final param = _in_path0(source, true, opts)
        inputs.addParam(param)
    }

    private String _in_path0(Object source, boolean pathQualifier, Map opts) {
        if( !opts.stageAs && opts.name )
            opts.stageAs = opts.remove('name')

        if( source instanceof LazyVar ) {
            final var = (LazyVar)source
            inputs.addFile(new ProcessFileInput(var, var.name, pathQualifier, opts))
            return var.name
        }
        else if( source instanceof CharSequence ) {
            final param = "\$in${inputs.size()}"
            if( !opts.stageAs )
                opts.stageAs = source
            inputs.addFile(new ProcessFileInput(new LazyVar(param), null, pathQualifier, opts))
            return param
        }
        else
            throw new IllegalArgumentException()
    }

    void _in_stdin() {
        final param = "\$in${inputs.size()}".toString()
        inputs.addParam(param)
        inputs.stdin = new LazyVar(param)
    }

    void _in_stdin(LazyVar var) {
        inputs.addParam(var.name)
        inputs.stdin = var
    }

    @CompileDynamic
    void _in_tuple(Object... elements) {
        if( elements.length < 2 )
            throw new IllegalArgumentException("Input `tuple` must define at least two elements -- Check process `$processName`")

        final param = "\$in${inputs.size()}".toString()
        inputs.addParam(param)

        for( int i = 0; i < elements.length; i++ ) {
            final item = elements[i]

            if( item instanceof LazyVar ) {
                final var = (LazyVar)item
                throw new IllegalArgumentException("Unqualified input value declaration is not allowed - replace `tuple ${var.name},..` with `tuple val(${var.name}),..`")
            }
            else if( item instanceof TokenValCall && item.val instanceof LazyVar ) {
                inputs.addVariable(item.val.name, new LazyTupleElement(param, i))
            }
            else if( item instanceof TokenEnvCall && item.val instanceof LazyVar ) {
                inputs.addEnv(item.val.name, new LazyTupleElement(param, i))
            }
            else if( item instanceof TokenFileCall ) {
                final name = _in_path0(item.target, false, [:])
                inputs.addVariable(name, new LazyTupleElement(param, i))
            }
            else if( item instanceof TokenPathCall ) {
                final name = _in_path0(item.target, true, item.opts)
                inputs.addVariable(name, new LazyTupleElement(param, i))
            }
            else if( item instanceof Map ) {
                throw new IllegalArgumentException("Unqualified input file declaration is not allowed - replace `tuple $item,..` with `tuple path(${item.key}, stageAs:'${item.value}'),..`")
            }
            else if( item instanceof GString ) {
                throw new IllegalArgumentException("Unqualified input file declaration is not allowed - replace `tuple \"$item\".. with `tuple path(\"$item\")..`")
            }
            else if( item instanceof TokenStdinCall || item == '-' ) {
                inputs.stdin = new LazyTupleElement(param, i)
            }
            else if( item instanceof String ) {
                throw new IllegalArgumentException("Unqualified input file declaration is not allowed - replace `tuple '$item',..` with `tuple path('$item'),..`")
            }
            else
                throw new IllegalArgumentException()
        }
    }

    void _in_val(LazyVar var) {
        inputs.addParam(var.name)
    }

    /// OUTPUTS

    void _out_env(Map opts=[:], Object target) {
        if( opts.emit )
            opts.name = opts.remove('emit')

        final name = _out_env0(target)
        outputs.addEnv(name)
        outputs.addParam(new LazyEnvCall(name), opts)
    }

    String _out_env0(Object target) {
        if( target instanceof LazyVar )
            return target.name
        else if( target instanceof CharSequence )
            return target.toString()
        else
            throw new IllegalArgumentException("Unexpected environment output definition - it should be either a string or a variable identifier - offending value: ${target?.getClass()?.getName()}")
    }

    void _out_eval(Map opts=[:], CharSequence cmd) {
        if( opts.emit )
            opts.name = opts.remove('emit')

        final name = outputs.addEval(cmd)
        outputs.addParam(new LazyEvalCall(name), opts)
    }

    void _out_file(Object target) {
        // note: check that is a String type to avoid to force
        // the evaluation of GString object to a string
        if( target instanceof String && target == '-' ) {
            _out_stdout()
            return
        }

        final key = _out_path0(target, false, [:])
        outputs.addParam(new LazyPathCall(key), [:])
    }

    void _out_path(Map opts=[:], Object target) {
        // note: check that is a String type to avoid to force
        // the evaluation of GString object to a string
        if( target instanceof String && target == '-' ) {
            _out_stdout(opts)
            return
        }

        // separate param options from path options
        final paramOpts = [optional: opts.optional]
        if( opts.emit )
            paramOpts.name = opts.remove('emit')

        final key = _out_path0(target, true, opts)
        outputs.addParam(new LazyPathCall(key), paramOpts)
    }

    private String _out_path0(Object target, boolean pathQualifier, Map opts) {
        outputs.addFile(new ProcessFileOutput(target, pathQualifier, opts))
    }

    void _out_stdout(Map opts=[:]) {
        if( opts.emit )
            opts.name = opts.remove('emit')

        outputs.addParam(new LazyVar('stdout'), opts)
    }

    @CompileDynamic
    void _out_tuple(Map opts=[:], Object... elements) {
        if( elements.length < 2 )
            throw new IllegalArgumentException("Output `tuple` must define at least two elements -- Check process `$processName`")

        // separate param options from path options
        if( opts.emit )
            opts.name = opts.remove('emit')

        // make lazy list with tuple elements
        final target = new LazyList(elements.size())

        for( int i = 0; i < elements.length; i++ ) {
            final item = elements[i]

            if( item instanceof LazyVar ) {
                throw new IllegalArgumentException("Unqualified output value declaration is not allowed - replace `tuple ${item.name},..` with `tuple val(${item.name}),..`")
            }
            else if( item instanceof TokenValCall ) {
                target << item.val
            }
            else if( item instanceof TokenEnvCall ) {
                final name = _out_env0(item.val)
                outputs.addEnv(name)
                target << new LazyEnvCall(name)
            }
            else if( item instanceof TokenEvalCall ) {
                final name = outputs.addEval(item.val)
                target << new LazyEvalCall(name)
            }
            else if( item instanceof TokenFileCall ) {
                // file pattern can be a String or GString
                final key = _out_path0(item.target, false, [optional: opts.optional])
                target << new LazyPathCall(key)
            }
            else if( item instanceof TokenPathCall ) {
                // file pattern can be a String or GString
                final key = _out_path0(item.target, true, item.opts + [optional: opts.optional])
                target << new LazyPathCall(key)
            }
            else if( item instanceof GString ) {
                throw new IllegalArgumentException("Unqualified output path declaration is not allowed - replace `tuple \"$item\",..` with `tuple path(\"$item\"),..`")
            }
            else if( item instanceof TokenStdoutCall || item == '-' ) {
                target << new LazyVar('stdout')
            }
            else if( item instanceof String ) {
                throw new IllegalArgumentException("Unqualified output path declaration is not allowed - replace `tuple '$item',..` with `tuple path('$item'),..`")
            }
            else
                throw new IllegalArgumentException("Invalid `tuple` output parameter declaration -- item: ${item}")
        }

        outputs.addParam(target, opts)
    }

    void _out_val(Map opts=[:], Object target) {
        outputs.addParam(target, opts)
    }

    /// BUILD

    ProcessDef build() {
        config.setInputs(inputs)
        config.setOutputs(outputs)
        super.build()
    }

}

@CompileStatic
class LazyTupleElement extends LazyVar {
    int index

    LazyTupleElement(String name, int index) {
        super(name)
        this.index = index
    }

    @Override
    Object resolve(Object binding) {
        final tuple = super.resolve(binding)
        if( tuple instanceof List )
            return tuple[index]
        else
            throw new IllegalArgumentException("Lazy binding of `${name}[${index}]` failed because `${name}` is not a tuple")
    }
}

@CompileStatic
class LazyEnvCall implements LazyAware {
    String name

    LazyEnvCall(String name) {
        this.name = name
    }

    @Override
    Object resolve(Object binding) {
        if( binding !instanceof TaskOutputCollector )
            throw new IllegalStateException()

        ((TaskOutputCollector)binding).env(name)
    }
}

@CompileStatic
class LazyEvalCall implements LazyAware {
    String name

    LazyEvalCall(String name) {
        this.name = name
    }

    @Override
    Object resolve(Object binding) {
        if( binding !instanceof TaskOutputCollector )
            throw new IllegalStateException()

        ((TaskOutputCollector)binding).eval(name)
    }
}

@CompileStatic
class LazyPathCall implements LazyAware {
    String key

    LazyPathCall(String key) {
        this.key = key
    }

    @Override
    Object resolve(Object binding) {
        if( binding !instanceof TaskOutputCollector )
            throw new IllegalStateException()

        ((TaskOutputCollector)binding).path(key)
    }
}
