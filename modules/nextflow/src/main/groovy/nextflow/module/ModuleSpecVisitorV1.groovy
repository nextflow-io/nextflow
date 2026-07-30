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
import nextflow.script.ast.ProcessNodeV1
import org.codehaus.groovy.ast.ClassHelper
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.NamedArgumentListExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.Statement

import static nextflow.module.ModuleSpec.ModuleParam
import static nextflow.script.ast.ASTUtils.*

/**
 * AST visitor to extract inputs/outputs from a legacy process.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ModuleSpecVisitorV1 {

    private ModuleSpec oldSpec

    ModuleSpecVisitorV1(ModuleSpec oldSpec) {
        this.oldSpec = oldSpec
    }

    List<ModuleParam> visitInputs(ProcessNodeV1 node) {
        return moduleInputs(asBlockStatements(node.inputs), oldSpec.inputs)
    }

    List<ModuleParam> visitOutputs(ProcessNodeV1 node) {
        return moduleOutputs(asBlockStatements(node.outputs), oldSpec.outputs)
    }

    List<ModuleParam> visitTopics(ProcessNodeV1 node) {
        return moduleTopics(asBlockStatements(node.outputs), oldSpec.topics)
    }

    private static List<ModuleParam> moduleInputs(List<Statement> statements, List<ModuleParam> oldParams) {
        final result = new ArrayList<ModuleParam>(statements.size())
        for( int i = 0; i < statements.size(); i++ ) {
            final call = asMethodCallX(statements[i])
            final oldParam = oldParams && i < oldParams.size() ? oldParams[i] : null
            result << moduleParam(call, null, oldParam)
        }
        return result
    }

    private static List<ModuleParam> moduleOutputs(List<Statement> statements, List<ModuleParam> oldParams) {
        final result = new ArrayList<ModuleParam>(statements.size())
        for( int i = 0; i < statements.size(); i++ ) {
            final call = asMethodCallX(statements[i])
            final emitName = outputName(call, 'emit')
            final topicName = outputName(call, 'topic')
            if( !emitName && topicName )
                continue
            final j = result.size()
            final oldParam = oldParams && j < oldParams.size() ? oldParams[j] : null
            result << moduleParam(call, emitName, oldParam)
        }
        return result
    }

    private static List<ModuleParam> moduleTopics(List<Statement> statements, List<ModuleParam> oldTopics) {
        final result = new ArrayList<ModuleParam>(statements.size())
        for( int i = 0; i < statements.size(); i++ ) {
            final call = asMethodCallX(statements[i])
            final topicName = outputName(call, 'topic')
            if( !topicName )
                continue
            final j = result.size()
            final oldParam = oldTopics && j < oldTopics.size() ? oldTopics[j] : null
            result << moduleParam(call, null, oldParam)
        }
        return result
    }

    private static String outputName(MethodCallExpression output, String type) {
        return Optional.of(output)
            .flatMap(call -> Optional.ofNullable(asNamedArgs(call)))
            .flatMap(namedArgs ->
                namedArgs.stream()
                    .filter(entry -> entry.getKeyExpression().getText() == type)
                    .findFirst()
            )
            .flatMap(entry -> Optional.ofNullable(
                entry.getValueExpression() instanceof VariableExpression
                    ? entry.getValueExpression().getName()
                    : null
            ))
            .orElse(null)
    }

    private static ModuleParam moduleParam(MethodCallExpression call, String emitName, ModuleParam oldParam) {
        final qualifier = call.getMethodAsString()
        final arguments = positionalArgs(call)
        if( qualifier == 'tuple' )
            return moduleTupleParam(arguments, oldParam)
        final paramName = paramName(arguments) ?: emitName
        final paramType = paramType(qualifier, arguments)
        return new ModuleParam(
            name: paramName ?: oldParam?.name,
            type: paramType ?: oldParam?.type,
            description: oldParam?.description,
            _passthrough: oldParam?._passthrough
        )
    }

    private static List<Expression> positionalArgs(MethodCallExpression call) {
        def arguments = asMethodCallArguments(call)
        return arguments[0] instanceof NamedArgumentListExpression
            ? arguments.tail()
            : arguments
    }

    private static ModuleParam moduleTupleParam(List<Expression> arguments, ModuleParam oldParam) {
        final components = new ArrayList<ModuleParam>(arguments.size())
        for( int i = 0; i < arguments.size(); i++ ) {
            final arg = (MethodCallExpression) arguments[i]
            final oldComponent = oldParam?.isTuple() && i < oldParam.components.size() ? oldParam.components[i] : null
            components << moduleParam(arg, null, oldComponent)
        }
        return new ModuleParam(
            name: null,
            type: null,
            components: components
        )
    }

    private static String paramName(List<Expression> arguments) {
        if( arguments.isEmpty() )
            return null
        final lastArg = arguments.last()
        return lastArg instanceof VariableExpression
            ? lastArg.getName()
            : null
    }

    private static String paramType(String qualifier, List<Expression> arguments) {
        switch( qualifier ) {
            case 'val':
                return valType(arguments.last())
            case 'file':
            case 'path':
                return 'file'
            case 'env':
            case 'eval':
            case 'stdin':
            case 'stdout':
                return 'string'
            default:
                return null
        }
    }

    private static String valType(Expression arg) {
        if( arg instanceof VariableExpression && arg.getName() == 'meta' )
            return 'map'
        switch( arg.getType() ) {
            case ClassHelper.Boolean_TYPE:
                return 'boolean'
            case ClassHelper.Float_TYPE:
                return 'float'
            case ClassHelper.Integer_TYPE:
                return 'integer'
            case ClassHelper.LIST_TYPE:
                return 'list'
            case ClassHelper.MAP_TYPE:
                return 'map'
            case ClassHelper.STRING_TYPE:
            case ClassHelper.GSTRING_TYPE:
                return 'string'
            default:
                return null
        }
    }

}
