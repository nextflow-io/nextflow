/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import java.lang.annotation.ElementType
import java.lang.annotation.Retention
import java.lang.annotation.RetentionPolicy
import java.lang.annotation.Target

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.control.CompilePhase
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.transform.ASTTransformation
import org.codehaus.groovy.transform.GroovyASTTransformation
import org.codehaus.groovy.transform.GroovyASTTransformationClass

/**
 * Declares an AST xform to escape and manipulate task command
 * special value
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Retention(RetentionPolicy.SOURCE)
@Target(ElementType.METHOD)
@GroovyASTTransformationClass(classes = [TaskCmdXformImpl])
@interface TaskCmdXform {

    //--- === implementation === ---

    @Slf4j
    @CompileStatic
    @GroovyASTTransformation(phase = CompilePhase.CONVERSION)
    class TaskCmdXformImpl implements ASTTransformation {
        @Override
        void visit(ASTNode[] nodes, SourceUnit source) {
            new TaskCmdXformVisitor(source).visitClass((ClassNode)nodes[1])
        }
    }
}
