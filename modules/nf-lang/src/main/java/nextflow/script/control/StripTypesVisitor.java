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

package nextflow.script.control;

import java.lang.reflect.Modifier;
import java.util.List;

import nextflow.script.ast.RecordNode;
import nextflow.script.types.Channel;
import nextflow.script.types.Record;
import nextflow.script.types.Tuple;
import nextflow.script.types.Value;
import org.codehaus.groovy.ast.ClassCodeExpressionTransformer;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.FieldNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.CastExpression;
import org.codehaus.groovy.ast.expr.ClassExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.DeclarationExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.syntax.Types;

import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Strip type annotations that are used by the Nextflow type checker
 * but not supported by the Groovy runtime:
 *
 * - The dataflow types (Channel and Value) are not compatible with
 *   the runtime types used when static typing is disabled.
 *
 * - The Tuple type can be specified with variable type arguments,
 *   but this is not supported by the JVM.
 *
 * - User-defined record types are only used to validate records with
 *   type RecordMap. Record types must be stripped from type annotations
 *   and type casts, and instanceof expressions on a record type must be
 *   replaced with a runtime check.
 *
 * @author Ben Sherman <bentshermman@gmail.com>
 */
public class StripTypesVisitor extends ClassCodeExpressionTransformer {

    private static final List<ClassNode> STRIP_TYPES = List.of(
        ClassHelper.makeWithoutCaching("nextflow.Channel"),
        ClassHelper.makeCached(Channel.class),
        ClassHelper.makeCached(Record.class),
        ClassHelper.makeCached(Tuple.class),
        ClassHelper.makeCached(Value.class)
    );

    private SourceUnit sourceUnit;

    public StripTypesVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    @Override
    public void visitMethod(MethodNode node) {
        // Erase record type parameters so that records (with type Record)
        // can be dispatched to these methods at runtime
        for( var param : node.getParameters() ) {
            if( isNamedRecordType(param.getType()) )
                param.setType(ClassHelper.dynamicType());
        }
        super.visitMethod(node);
    }

    @Override
    public Expression transform(Expression node) {
        if( node instanceof BinaryExpression be && be.getOperation().getType() == Types.KEYWORD_INSTANCEOF ) {
            return transformInstanceof(be);
        }

        if( node instanceof CastExpression ce ) {
            return transformCast(ce);
        }

        if( node instanceof ClosureExpression ce ) {
            ce.visit(this);
            return ce;
        }

        if( node instanceof DeclarationExpression de ) {
            stripTypeAnnotation(de);
        }

        return super.transform(node);
    }

    private Expression transformCast(CastExpression node) {
        var type = node.getType();
        if( type.getGenericsTypes() != null ) {
            var fn = parameterizedType(type);
            return callThisX("_as_type", args(transform(node.getExpression()), classX(fn.getDeclaringClass()), constX(fn.getName())));
        }
        else {
            return callThisX("_as_type", args(transform(node.getExpression()), classX(type)));
        }
    }

    private ClassNode hiddenClassNode;

    private FieldNode parameterizedType(ClassNode type) {
        if( hiddenClassNode == null ) {
            var moduleNode = sourceUnit.getAST();
            var scriptClass = moduleNode.getClasses().get(0);
            var packageName = scriptClass.getNameWithoutPackage();
            this.hiddenClassNode = new RecordNode(packageName + "." + "__ParameterizedTypes");
            moduleNode.addClass(hiddenClassNode);
        }
        var name = "_" + hiddenClassNode.getFields().size();
        var fn = new FieldNode(
            name,
            Modifier.PUBLIC,
            type,
            hiddenClassNode,
            null
        );
        hiddenClassNode.addField(fn);
        return fn;
    }

    private Expression stripTypeAnnotation(DeclarationExpression node) {
        if( node.getLeftExpression() instanceof VariableExpression ve ) {
            if( shouldStripType(ve.getType()) )
                node.setLeftExpression(new VariableExpression(ve.getName()));
        }
        return node;
    }

    private Expression transformInstanceof(BinaryExpression node) {
        var right = node.getRightExpression();
        if( right instanceof ClassExpression ce && isNamedRecordType(ce.getType()) )
            return callThisX("_instanceof_record_type", args(transform(node.getLeftExpression()), right));
        return super.transform(node);
    }

    private boolean shouldStripType(ClassNode type) {
        return STRIP_TYPES.contains(type) || isNamedRecordType(type);
    }

    private boolean isNamedRecordType(ClassNode type) {
        return type.redirect() instanceof RecordNode;
    }

}
