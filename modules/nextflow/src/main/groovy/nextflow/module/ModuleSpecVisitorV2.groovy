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

package nextflow.module

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.script.ast.AssignmentExpression
import nextflow.script.ast.ProcessNodeV2
import nextflow.script.ast.TupleParameter
import org.codehaus.groovy.ast.ClassHelper
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.Parameter
import org.codehaus.groovy.ast.expr.BinaryExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.Statement

import static nextflow.module.ModuleSpec.ModuleParam
import static nextflow.script.ast.ASTUtils.*

/**
 * AST visitor to extract inputs/outputs from a typed process.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ModuleSpecVisitorV2 {

    private ModuleSpec oldSpec

    ModuleSpecVisitorV2(ModuleSpec oldSpec) {
        this.oldSpec = oldSpec
    }

    List<ModuleParam> visitInputs(ProcessNodeV2 node) {
        return moduleInputs(node.inputs, oldSpec.inputs)
    }

    List<ModuleParam> visitOutputs(ProcessNodeV2 node) {
        return moduleOutputs(asBlockStatements(node.outputs), oldSpec.outputs)
    }

    List<ModuleParam> visitTopics(ProcessNodeV2 node) {
        return moduleTopics(asBlockStatements(node.topics), oldSpec.topics)
    }

    private static List<ModuleParam> moduleInputs(Parameter[] params, List<ModuleParam> oldParams) {
        final result = new ArrayList<ModuleParam>(params.length)
        for( int i = 0; i < params.length; i++ ) {
            final param = params[i]
            final oldParam = oldParams && i < oldParams.size() ? oldParams[i] : null
            if( param instanceof TupleParameter ) {
                result << moduleTupleInput(param, oldParam)
            } else {
                result << moduleInput(param, oldParam)
            }
        }
        return result
    }

    private static ModuleParam moduleTupleInput(TupleParameter tp, ModuleParam oldParam) {
        final components = new ArrayList<ModuleParam>(tp.components.length)
        for( int i = 0; i < tp.components.length; i++ ) {
            final el = tp.components[i]
            final oldComponent = oldParam?.isTuple() && i < oldParam.components.size() ? oldParam.components[i] : null
            components << moduleInput(el, oldComponent)
        }
        return new ModuleParam(
            name: null,
            type: null,
            components: components
        )
    }

    private static ModuleParam moduleInput(Parameter param, ModuleParam oldParam) {
        return new ModuleParam(
            name: param.getName(),
            type: paramType(param.getType()) ?: oldParam?.type,
            description: oldParam?.description,
            _passthrough: oldParam?._passthrough
        )
    }

    private static List<ModuleParam> moduleOutputs(List<Statement> statements, List<ModuleParam> oldParams) {
        final result = new ArrayList<ModuleParam>(statements.size())

        for( int i = 0; i < statements.size(); i++ ) {
            final output = ((ExpressionStatement) statements[i]).getExpression()
            final oldParam = oldParams && i < oldParams.size() ? oldParams[i] : null
            if( isRecordOutput(output) )
                result << moduleRecordOutput((MethodCallExpression) output, oldParam)
            else if( isTupleOutput(output) )
                result << moduleTupleOutput((MethodCallExpression) output, oldParam)
            else
                result << moduleOutput(output, oldParam)
        }

        return result
    }

    private static boolean isRecordOutput(Expression output) {
        return output instanceof MethodCallExpression
            && output.getMethodAsString() == 'record'
    }

    private static ModuleParam moduleRecordOutput(MethodCallExpression output, ModuleParam oldParam) {
        final namedArgs = asNamedArgs(output)
        final components = new ArrayList<ModuleParam>(namedArgs.size())
        for( int i = 0; i < namedArgs.size(); i++ ) {
            final entry = namedArgs[i]
            final oldComponent = oldParam?.isTuple() && i < oldParam.components.size() ? oldParam.components[i] : null
            final name = entry.getKeyExpression().getText()
            final type = paramType(entry.getValueExpression())

            final component = new ModuleParam(
                name: name ?: oldParam?.name,
                type: type ?: oldParam?.type,
                description: oldParam?.description,
                _passthrough: oldParam?._passthrough
            )
            components << component
        }
        return new ModuleParam(
            name: null,
            type: null,
            components: components
        )
    }

    private static boolean isTupleOutput(Expression output) {
        return output instanceof MethodCallExpression
            && output.getMethodAsString() == 'tuple'
    }

    private static ModuleParam moduleTupleOutput(MethodCallExpression output, ModuleParam oldParam) {
        final arguments = asMethodCallArguments(output)
        final components = new ArrayList<ModuleParam>(arguments.size())
        for( int i = 0; i < arguments.size(); i++ ) {
            final el = arguments[i]
            final oldComponent = oldParam?.isTuple() && i < oldParam.components.size() ? oldParam.components[i] : null
            components << moduleOutput(el, oldComponent)
        }
        return new ModuleParam(
            name: null,
            type: null,
            components: components
        )
    }

    private static ModuleParam moduleOutput(Expression output, ModuleParam oldParam) {
        final target = 
            output instanceof AssignmentExpression ? (VariableExpression) output.getLeftExpression() : 
            output instanceof VariableExpression ? output : null

        final name = target != null ? target.getName() : null
        final type = target != null ? paramType(target.getType()) : paramType(output)

        return new ModuleParam(
            name: name ?: oldParam?.name,
            type: type ?: oldParam?.type,
            description: oldParam?.description,
            _passthrough: oldParam?._passthrough
        )
    }

    private static List<ModuleParam> moduleTopics(List<Statement> statements, List<ModuleParam> oldTopics) {
        final result = new ArrayList<ModuleParam>(statements.size())

        for( int i = 0; i < statements.size(); i++ ) {
            final expr = ((ExpressionStatement) statements[i]).getExpression()
            final be = (BinaryExpression) expr
            final source = be.getLeftExpression()
            final oldParam = oldTopics && i < oldTopics.size() ? oldTopics[i] : null
            if( isRecordOutput(source) )
                result << moduleRecordOutput((MethodCallExpression) source, oldParam)
            else if( isTupleOutput(source) )
                result << moduleTupleOutput((MethodCallExpression) source, oldParam)
            else
                result << moduleOutput(source, oldParam)
        }

        return result
    }

    private static final ClassNode PATH_TYPE = ClassHelper.makeCached(java.nio.file.Path)

    private static String paramType(ClassNode type) {
        if( !type || !type.isResolved() )
            return null

        if( type.implementsInterface(ClassHelper.ITERABLE_TYPE) && !type.equals(PATH_TYPE) ) {
            final gts = type.getGenericsTypes()
            final elementType = gts ? gts[0].getType() : null
            return PATH_TYPE.equals(elementType) ? 'file' : 'list'
        }

        switch( type ) {
            case ClassHelper.Boolean_TYPE:
                return 'boolean'
            case PATH_TYPE:
                return 'file'
            case ClassHelper.Float_TYPE:
                return 'float'
            case ClassHelper.Integer_TYPE:
                return 'integer'
            case ClassHelper.MAP_TYPE:
                return 'map'
            case ClassHelper.STRING_TYPE:
            case ClassHelper.GSTRING_TYPE:
                return 'string'
            default:
                return null
        }
    }

    private static String paramType(Expression node) {
        if( node instanceof MethodCallExpression ) {
            final name = node.getMethodAsString()
            switch( name ) {
                case 'env':
                case 'eval':
                case 'stdout':
                    return 'string'
                case 'file':
                case 'files':
                    return 'file'
            }
        }

        return paramType(node.getType())
    }

}
