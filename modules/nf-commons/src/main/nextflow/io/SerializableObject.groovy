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

package nextflow.io

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
 * Declares an AST xform to automatically add the {@link SerializableMarker} interface
 * for a class marked as @SerializableObject
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Retention(RetentionPolicy.RUNTIME)
@Target([ElementType.TYPE])
@GroovyASTTransformationClass(classes = [ValueObjectImpl])
@interface SerializableObject {

    //--- === implementation === ---

    @Slf4j
    @CompileStatic
    @GroovyASTTransformation(phase = CompilePhase.SEMANTIC_ANALYSIS)
    class ValueObjectImpl implements ASTTransformation {
        @Override
        void visit(ASTNode[] nodes, SourceUnit source) {
            final target = (ClassNode)nodes[1]
            target.addInterface( new ClassNode(SerializableMarker.class) )
        }
    }
}
