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
package nextflow.antlr

import groovy.transform.CompileStatic
import groovy.transform.ToString
import org.codehaus.groovy.ast.MethodNode
import org.codehaus.groovy.ast.stmt.Statement

/**
 * Container for Nextflow AST elements which will be
 * injected into the Groovy script class.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@ToString
class Module {
    private List<MethodNode> methods = []
    private List<Statement> statements = []

    void addMethod(MethodNode method) {
        methods.add(method)
    }

    List<MethodNode> getMethods() {
        methods
    }

    void addStatement(Statement stmt) {
        statements.add(stmt)
    }

    List<Statement> getStatements() {
        statements
    }
}
