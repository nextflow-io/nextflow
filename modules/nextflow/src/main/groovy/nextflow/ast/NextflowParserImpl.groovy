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

package nextflow.ast

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.script.BaseScript
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.MethodNode
import org.codehaus.groovy.control.CompilePhase
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.syntax.SyntaxException
import org.codehaus.groovy.transform.ASTTransformation
import org.codehaus.groovy.transform.GroovyASTTransformation

@Slf4j
@CompileStatic
@GroovyASTTransformation(phase = CompilePhase.CONVERSION)
class NextflowParserImpl implements ASTTransformation {

    static private Set<String> RESERVED_NAMES

    static {
        // method names implicitly defined by the groovy script SHELL
        RESERVED_NAMES = ['main','run','runScript'] as Set
        // existing method cannot be used for custom script definition
        for( def method : BaseScript.getMethods() ) {
            RESERVED_NAMES.add(method.name)
        }
    }

    @Override
    void visit(ASTNode[] astNodes, SourceUnit unit) {
        new ParserCodeVisitor(unit).visitClass((ClassNode)astNodes[1])
    }

    static class ParserCodeVisitor extends ClassCodeVisitorSupport {

        private SourceUnit unit

        private Set<String> processNames = []

        private Set<String> workflowNames = []

        private Set<String> functionNames = []

        ParserCodeVisitor(SourceUnit unit) {
            this.unit = unit
        }

        @Override
        protected SourceUnit getSourceUnit() { unit }

        @Override
        void visitMethod(MethodNode method) {
            if( method.public && !method.static && !method.synthetic && !method.metaDataMap?.'org.codehaus.groovy.ast.MethodNode.isScriptBody') {
                if( !isIllegalName(method.name, method))
                    functionNames.add(method.name)
            }

            // TODO: fix lazy gstring
            // TODO: fix output emit/topic options
            // TODO: wrap output property expression in closure
        }

        protected boolean isIllegalName(String name, ASTNode node) {
            if( name in RESERVED_NAMES ) {
                unit.addError( new SyntaxException("Identifier `$name` is reserved for internal use", node.lineNumber, node.columnNumber+8) )
                return true
            }
            if( name in workflowNames || name in processNames ) {
                unit.addError( new SyntaxException("Identifier `$name` is already used by another definition", node.lineNumber, node.columnNumber+8) )
                return true
            }
            if( name.contains(Const.SCOPE_SEP) ) {
                final offset = 8 + 2 + name.indexOf(Const.SCOPE_SEP)
                unit.addError( new SyntaxException("Process and workflow names cannot contain colon character", node.lineNumber, node.columnNumber+offset) )
                return true
            }
            return false
        }

    }

}
