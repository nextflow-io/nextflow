/*
 * Copyright 2024-2025, Seqera Labs
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

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import groovy.lang.Tuple2;
import nextflow.script.ast.ASTNodeMarker;
import org.codehaus.groovy.GroovyBugError;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassCodeExpressionTransformer;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.DynamicVariable;
import org.codehaus.groovy.ast.GenericsType;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.ClassExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.ConstructorCallExpression;
import org.codehaus.groovy.ast.expr.DeclarationExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.CatchStatement;
import org.codehaus.groovy.control.ClassNodeResolver;
import org.codehaus.groovy.control.CompilationUnit;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.runtime.memoize.UnlimitedConcurrentCache;
import org.codehaus.groovy.syntax.SyntaxException;
import org.codehaus.groovy.syntax.Types;
import org.codehaus.groovy.vmplugin.VMPluginFactory;

import static groovy.lang.Tuple.tuple;
import static org.codehaus.groovy.ast.tools.ClosureUtils.getParametersSafe;

/**
 * Resolve the names of symbols.
 *
 * @see org.codehaus.groovy.control.ResolveVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ResolveVisitor extends ClassCodeExpressionTransformer {

    public static final String[] DEFAULT_PACKAGE_PREFIXES = { "java.lang.", "java.util.", "java.io.", "java.net.", "groovy.lang.", "groovy.util." };

    public static final String[] EMPTY_STRING_ARRAY = new String[0];

    private SourceUnit sourceUnit;

    private CompilationUnit compilationUnit;

    private ClassNodeResolver classNodeResolver = new ClassNodeResolver();

    /**
     * Default imports can be accessed only by their simple name.
     */
    private List<ClassNode> defaultImports;

    /**
     * Lib imports can be accessed only by their fully-qualified name.
     */
    private List<ClassNode> libImports;

    public ResolveVisitor(SourceUnit sourceUnit, CompilationUnit compilationUnit, List<ClassNode> defaultImports, List<ClassNode> libImports) {
        this.sourceUnit = sourceUnit;
        this.compilationUnit = compilationUnit;
        this.defaultImports = defaultImports;
        this.libImports = libImports;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    @Override
    public void visitCatchStatement(CatchStatement cs) {
        resolveOrFail(cs.getExceptionType(), cs);
        if( ClassHelper.isDynamicTyped(cs.getExceptionType()) )
            cs.getVariable().setType(ClassHelper.make(Exception.class));
        super.visitCatchStatement(cs);
    }

    public void resolveOrFail(ClassNode type, ASTNode node) {
        if( !resolve(type) )
            addError("`" + type.toString(false) + "` is not defined", node);
    }

    public boolean resolve(ClassNode type) {
        var genericsTypes = type.getGenericsTypes();
        resolveGenericsTypes(genericsTypes);

        if( type.isPrimaryClassNode() )
            return true;
        if( type.isResolved() )
            return true;
        if( resolveFromModule(type) )
            return true;
        if( resolveFromLibImports(type) )
            return true;
        if( !type.hasPackageName() && resolveFromDefaultImports(type) )
            return true;
        return resolveFromClassResolver(type.getName()) != null;
    }

    private boolean resolveGenericsTypes(GenericsType[] types) {
        if( types == null )
            return true;
        boolean resolved = true;
        for( var type : types ) {
            if( !resolveGenericsType(type) )
                resolved = false;
        }
        return resolved;
    }

    private boolean resolveGenericsType(GenericsType genericsType) {
        if( genericsType.isResolved() )
            return true;
        var type = genericsType.getType();
        resolveOrFail(type, genericsType);
        if( resolveGenericsTypes(type.getGenericsTypes()) )
            genericsType.setResolved(genericsType.getType().isResolved());
        return genericsType.isResolved();
    }

    protected boolean resolveFromModule(ClassNode type) {
        var name = type.getName();
        var module = sourceUnit.getAST();
        for( var cn : module.getClasses() ) {
            if( name.equals(cn.getName()) ) {
                if( cn != type )
                    type.setRedirect(cn);
                return true;
            }
        }
        return false;
    }

    protected boolean resolveFromLibImports(ClassNode type) {
        var name = type.getName();
        for( var cn : libImports ) {
            if( name.equals(cn.getName()) ) {
                type.setRedirect(cn);
                return true;
            }
        }
        return false;
    }

    protected boolean resolveFromDefaultImports(ClassNode type) {
        // resolve from script imports
        var typeName = type.getName();
        for( var cn : defaultImports ) {
            if( typeName.equals(cn.getNameWithoutPackage()) ) {
                type.setRedirect(cn);
                return true;
            }
        }
        // resolve from default imports cache
        var packagePrefixSet = DEFAULT_IMPORT_CLASS_AND_PACKAGES_CACHE.get(typeName);
        if( packagePrefixSet != null ) {
            if( resolveFromDefaultImports(type, packagePrefixSet.toArray(EMPTY_STRING_ARRAY)) )
                return true;
        }
        // resolve from default imports
        if( resolveFromDefaultImports(type, DEFAULT_PACKAGE_PREFIXES) ) {
            return true;
        }
        if( "BigInteger".equals(typeName) ) {
            type.setRedirect(ClassHelper.BigInteger_TYPE);
            return true;
        }
        if( "BigDecimal".equals(typeName) ) {
            type.setRedirect(ClassHelper.BigDecimal_TYPE);
            return true;
        }
        return false;
    }

    private static final Map<String, Set<String>> DEFAULT_IMPORT_CLASS_AND_PACKAGES_CACHE = new UnlimitedConcurrentCache<>();
    static {
        DEFAULT_IMPORT_CLASS_AND_PACKAGES_CACHE.putAll(VMPluginFactory.getPlugin().getDefaultImportClasses(DEFAULT_PACKAGE_PREFIXES));
    }

    protected boolean resolveFromDefaultImports(ClassNode type, String[] packagePrefixes) {
        var typeName = type.getName();
        for( var packagePrefix : packagePrefixes ) {
            var redirect = resolveFromClassResolver(packagePrefix + typeName);
            if( redirect != null ) {
                type.setRedirect(redirect);
                // don't update cache when using a cached lookup
                if( packagePrefixes == DEFAULT_PACKAGE_PREFIXES ) {
                    var packagePrefixSet = DEFAULT_IMPORT_CLASS_AND_PACKAGES_CACHE.computeIfAbsent(typeName, key -> new HashSet<>(2));
                    packagePrefixSet.add(packagePrefix);
                }
                return true;
            }
        }
        return false;
    }

    protected ClassNode resolveFromClassResolver(String name) {
        var lookupResult = classNodeResolver.resolveName(name, compilationUnit);
        if( lookupResult == null )
            return null;
        if( !lookupResult.isClassNode() )
            throw new GroovyBugError("class resolver lookup result is not a class node");
        return lookupResult.getClassNode();
    }

    @Override
    public Expression transform(Expression exp) {
        if( exp == null )
            return null;
        Expression result;
        if( exp instanceof VariableExpression ve ) {
            result = transformVariableExpression(ve);
        }
        else if( exp instanceof PropertyExpression pe ) {
            result = transformPropertyExpression(pe);
        }
        else if( exp instanceof DeclarationExpression de ) {
            result = transformDeclarationExpression(de);
        }
        else if( exp instanceof BinaryExpression be ) {
            result = transformBinaryExpression(be);
        }
        else if( exp instanceof MethodCallExpression mce ) {
            result = transformMethodCallExpression(mce);
        }
        else if( exp instanceof ClosureExpression ce ) {
            result = transformClosureExpression(ce);
        }
        else if( exp instanceof ConstructorCallExpression cce ) {
            result = transformConstructorCallExpression(cce);
        }
        else {
            resolveOrFail(exp.getType(), exp);
            result = exp.transformExpression(this);
        }
        if( result != null && result != exp ) {
            result.setSourcePosition(exp);
        }
        return result;
    }

    protected Expression transformVariableExpression(VariableExpression ve) {
        var v = ve.getAccessedVariable();
        if( v instanceof DynamicVariable ) {
            // attempt to resolve variable as type name
            var name = ve.getName();
            var type = ClassHelper.make(name);
            var isClass = type.isResolved();
            if( !isClass )
                isClass = resolve(type);
            if( isClass )
                return new ClassExpression(type);
        }
        if( inVariableDeclaration ) {
            // resolve type of variable declaration
            resolveOrFail(ve.getType(), ve);
            var origin = ve.getOriginType();
            if( origin != ve.getType() )
                resolveOrFail(origin, ve);
        }
        // if the variable is still dynamic (i.e. unresolved), it will be handled by DynamicVariablesVisitor
        return ve;
    }

    protected Expression transformPropertyExpression(PropertyExpression pe) {
        Expression objectExpression;
        Expression property;
        objectExpression = transform(pe.getObjectExpression());
        property = transform(pe.getProperty());
        var result = new PropertyExpression(objectExpression, property, pe.isSafe());
        result.setSpreadSafe(pe.isSpreadSafe());
        // attempt to resolve property expression as a fully-qualified class name
        var className = lookupClassName(result);
        if( className != null ) {
            var type = ClassHelper.make(className);
            type.putNodeMetaData(ASTNodeMarker.FULLY_QUALIFIED, true);
            if( resolve(type) )
                return new ClassExpression(type);
        }
        return result;
    }

    private static String lookupClassName(PropertyExpression node) {
        boolean doInitialClassTest = true;
        StringBuilder name = new StringBuilder(32);
        Expression expr = node;
        while( expr != null && name != null ) {
            if( expr instanceof VariableExpression ve ) {
                var varName = ve.getName();
                var classNameInfo = makeClassName(doInitialClassTest, name, varName);
                name = classNameInfo.getV1();
                doInitialClassTest = classNameInfo.getV2();
                break;
            }

            if( expr instanceof PropertyExpression pe ) {
                var property = pe.getPropertyAsString();
                var classNameInfo = makeClassName(doInitialClassTest, name, property);
                name = classNameInfo.getV1();
                doInitialClassTest = classNameInfo.getV2();
                expr = pe.getObjectExpression();
            }
            else {
                return null;
            }
        }
        if( name == null || name.length() == 0 )
            return null;
        return name.toString();
    }

    private static Tuple2<StringBuilder, Boolean> makeClassName(boolean doInitialClassTest, StringBuilder name, String varName) {
        if( doInitialClassTest ) {
            return isValidClassName(varName)
                ? tuple(new StringBuilder(varName), Boolean.FALSE)
                : tuple(null, Boolean.TRUE);
        }
        name.insert(0, varName + ".");
        return tuple(name, Boolean.FALSE);
    }

    private static boolean isValidClassName(String name) {
        if( name == null || name.length() == 0 )
            return false;
        return !Character.isLowerCase(name.charAt(0));
    }

    private boolean inVariableDeclaration;

    protected Expression transformDeclarationExpression(DeclarationExpression de) {
        inVariableDeclaration = true;
        var left = transform(de.getLeftExpression());
        inVariableDeclaration = false;
        if( left instanceof ClassExpression ) {
            addError("`" + left.getType().getName() + "` is already defined as a type", de.getLeftExpression());
            return de;
        }
        var right = transform(de.getRightExpression());
        if( right == de.getRightExpression() )
            return de;
        var result = new DeclarationExpression(left, de.getOperation(), right);
        result.setDeclaringClass(de.getDeclaringClass());
        return result;
    }

    protected Expression transformBinaryExpression(BinaryExpression be) {
        var left = transform(be.getLeftExpression());
        if( Types.isAssignment(be.getOperation().getType()) && left instanceof ClassExpression ) {
            addError("`" + left.getType().getName() + "` is already defined as a type", be.getLeftExpression());
            return be;
        }
        be.setLeftExpression(left);
        be.setRightExpression(transform(be.getRightExpression()));
        return be;
    }

    protected Expression transformMethodCallExpression(MethodCallExpression mce) {
        var args = transform(mce.getArguments());
        var method = transform(mce.getMethod());
        var object = mce.getObjectExpression();
        if( !mce.isImplicitThis() )
            object = transform(object);
        var result = new MethodCallExpression(object, method, args);
        result.setMethodTarget(mce.getMethodTarget());
        result.setImplicitThis(mce.isImplicitThis());
        result.setSpreadSafe(mce.isSpreadSafe());
        result.setSafe(mce.isSafe());
        result.copyNodeMetaData(mce);
        return result;
    }

    protected Expression transformClosureExpression(ClosureExpression ce) {
        for( var param : getParametersSafe(ce) ) {
            resolveOrFail(param.getType(), ce);
            if( param.hasInitialExpression() )
                param.setInitialExpression(transform(param.getInitialExpression()));
        }
        visit(ce.getCode());
        return ce;
    }

    protected Expression transformConstructorCallExpression(ConstructorCallExpression cce) {
        var cceType = cce.getType();
        resolveOrFail(cceType, cce);
        if( cceType.isAbstract() )
            addError("`" + cceType.getName() + "` is an abstract type and cannot be constructed directly", cce);
        return cce.transformExpression(this);
    }

    @Override
    public void addError(String message, ASTNode node) {
        var cause = new UnresolvedNameError(message, node);
        var errorMessage = new SyntaxErrorMessage(cause, sourceUnit);
        sourceUnit.getErrorCollector().addErrorAndContinue(errorMessage);
    }

    private class UnresolvedNameError extends SyntaxException implements PhaseAware {

        public UnresolvedNameError(String message, ASTNode node) {
            super(message, node);
        }

        @Override
        public int getPhase() {
            return Phases.NAME_RESOLUTION;
        }
    }

}
