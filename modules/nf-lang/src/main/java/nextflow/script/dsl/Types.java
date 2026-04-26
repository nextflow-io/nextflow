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
package nextflow.script.dsl;

import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.util.List;
import java.util.Map;

import nextflow.script.ast.ASTNodeMarker;
import org.codehaus.groovy.GroovyBugError;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.GenericsType;
import org.codehaus.groovy.ast.tools.GenericsUtils;

/**
 * Utility constants and functions for working with Nextflow types.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class Types {

    public static final List<ClassNode> DEFAULT_SCRIPT_IMPORTS = List.of(
        ClassHelper.makeCached(nextflow.script.types.Bag.class),
        ClassHelper.makeCached(nextflow.script.types.Channel.class),
        ClassHelper.makeCached(nextflow.script.types.Duration.class),
        ClassHelper.makeCached(nextflow.script.types.MemoryUnit.class),
        ClassHelper.makeCached(java.nio.file.Path.class),
        ClassHelper.makeCached(nextflow.script.types.Value.class),
        ClassHelper.makeCached(nextflow.script.types.VersionNumber.class)
    );

    public static final List<ClassNode> DEFAULT_CONFIG_IMPORTS = List.of(
        ClassHelper.makeCached(nextflow.script.types.Bag.class),
        ClassHelper.makeCached(nextflow.script.types.Duration.class),
        ClassHelper.makeCached(nextflow.script.types.MemoryUnit.class)
    );

    /**
     * Determine whether a type is a functional interface.
     *
     * @param type
     */
    public static boolean isFunctionalInterface(ClassNode type) {
        return type.getAnnotations().stream()
            .filter(an -> an.getClassNode().getName().equals(FunctionalInterface.class.getName()))
            .findFirst()
            .isPresent();
    }

    /**
     * Determine whether a class is a namespace.
     *
     * @param cn
     */
    public static boolean isNamespace(ClassNode cn) {
        return cn.implementsInterface(ClassHelper.makeCached(Namespace.class));
    }

    /**
     * Get the display name of a type.
     *
     * @param type
     */
    public static String getName(ClassNode type) {
        if( type == null )
            return "?";

        if( type.isArray() )
            return getName(type.getComponentType());

        if( isFunctionalInterface(type) )
            return closureName(type);

        return typeName(type);
    }

    private static String closureName(ClassNode type) {
        var mn = ClassHelper.findSAM(type);
        var spec = GenericsUtils.extractPlaceholders(type);
        var sb = new StringBuilder();

        var params = mn.getParameters();
        sb.append('(');
        for( int i = 0; i < params.length; i++ ) {
            if( i > 0 )
                sb.append(", ");
            var paramType = specificType(params[i].getType(), spec);
            sb.append(getName(paramType));
        }
        sb.append(')');

        var returnType = specificType(mn.getReturnType(), spec);
        sb.append(" -> ");
        sb.append(
            ClassHelper.VOID_TYPE.equals(returnType)
                ? "()"
                : getName(returnType)
        );

        return sb.toString();
    }

    private static ClassNode specificType(ClassNode type, Map<GenericsType.GenericsTypeName, GenericsType> spec) {
        if( !type.isGenericsPlaceHolder() )
            return type;
        var name = type.getUnresolvedName();
        return spec.get(new GenericsType.GenericsTypeName(name)).getType();
    }

    private static String tupleName(ClassNode type) {
        var sb = new StringBuilder();
        sb.append('(');
        genericsTypeNames(type.getGenericsTypes(), sb);
        sb.append(')');
        return sb.toString();
    }

    private static String typeName(ClassNode type) {
        var sb = new StringBuilder();

        var placeholder = type.isGenericsPlaceHolder();
        if( placeholder )
            sb.append(type.getUnresolvedName());
        else if( type.getNodeMetaData(ASTNodeMarker.FULLY_QUALIFIED) != null )
            sb.append(type.getName());
        else if( type.isResolved() )
            sb.append(getName(type.getTypeClass()));
        else
            sb.append(getName(type.getNameWithoutPackage()));

        if( !placeholder && type.getGenericsTypes() != null ) {
            sb.append('<');
            genericsTypeNames(type.getGenericsTypes(), sb);
            sb.append('>');
        }

        if( type.getNodeMetaData(ASTNodeMarker.NULLABLE) != null )
            sb.append('?');

        return sb.toString();
    }

    private static void genericsTypeNames(GenericsType[] genericsTypes, StringBuilder sb) {
        for( int i = 0; i < genericsTypes.length; i++ ) {
            if( i > 0 )
                sb.append(", ");
            sb.append(getName(genericsTypes[i].getType()));
        }
    }

    public static String getName(Type type) {
        return
            type instanceof Class c ? getName(c) :
            type instanceof ParameterizedType pt ? getName(pt) :
            getName(type.getTypeName());
    }

    private static String getName(Class type) {
        return getName(type.getSimpleName());
    }

    private static String getName(ParameterizedType pt) {
        var sb = new StringBuilder();
        sb.append(getName(pt.getRawType()));
        sb.append('<');
        for( int i = 0; i < pt.getActualTypeArguments().length; i++ ) {
            if( i > 0 )
                sb.append(", ");
            sb.append(getName(pt.getActualTypeArguments()[i]));
        }
        sb.append('>');
        return sb.toString();
    }

    public static String getName(String name) {
        if( "Object".equals(name) )
            return "?";
        return name;
    }

}
