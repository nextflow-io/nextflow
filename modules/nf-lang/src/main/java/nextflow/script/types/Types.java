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
package nextflow.script.types;

import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.dsl.Namespace;
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
        ClassHelper.makeCached(Bag.class),
        ClassHelper.makeCached(Channel.class),
        ClassHelper.makeCached(Duration.class),
        ClassHelper.makeCached(MemoryUnit.class),
        ClassHelper.makeCached(Path.class),
        ClassHelper.makeCached(Value.class),
        ClassHelper.makeCached(VersionNumber.class)
    );

    public static final List<ClassNode> DEFAULT_CONFIG_IMPORTS = List.of(
        ClassHelper.makeCached(Bag.class),
        ClassHelper.makeCached(Duration.class),
        ClassHelper.makeCached(MemoryUnit.class)
    );

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

        if( ClassHelper.isFunctionalInterface(type) )
            return closureName(type);

        return typeName(type);
    }

    private static String closureName(ClassNode type) {
        var mn = ClassHelper.findSAM(type);
        var spec = GenericsUtils.extractPlaceholders(type);
        var builder = new StringBuilder();

        var params = mn.getParameters();
        builder.append('(');
        for( int i = 0; i < params.length; i++ ) {
            if( i > 0 )
                builder.append(", ");
            var paramType = specificType(params[i].getType(), spec);
            builder.append(getName(paramType));
        }
        builder.append(')');

        var returnType = specificType(mn.getReturnType(), spec);
        builder.append(" -> ");
        builder.append(
            ClassHelper.VOID_TYPE.equals(returnType)
                ? "()"
                : getName(returnType)
        );

        return builder.toString();
    }

    private static ClassNode specificType(ClassNode type, Map<GenericsType.GenericsTypeName, GenericsType> spec) {
        if( !type.isGenericsPlaceHolder() )
            return type;
        var name = type.getUnresolvedName();
        return spec.get(new GenericsType.GenericsTypeName(name)).getType();
    }

    private static String tupleName(ClassNode type) {
        var builder = new StringBuilder();
        builder.append('(');
        genericsTypeNames(type.getGenericsTypes(), builder);
        builder.append(')');
        return builder.toString();
    }

    private static String typeName(ClassNode type) {
        var builder = new StringBuilder();

        var placeholder = type.isGenericsPlaceHolder();
        if( placeholder )
            builder.append(type.getUnresolvedName());
        else if( type.getNodeMetaData(ASTNodeMarker.FULLY_QUALIFIED) != null )
            builder.append(type.getName());
        else if( type.isResolved() )
            builder.append(getName(type.getTypeClass()));
        else
            builder.append(getName(type.getNameWithoutPackage()));

        if( !placeholder && type.getGenericsTypes() != null ) {
            builder.append('<');
            genericsTypeNames(type.getGenericsTypes(), builder);
            builder.append('>');
        }

        if( type.getNodeMetaData(ASTNodeMarker.NULLABLE) != null )
            builder.append('?');

        return builder.toString();
    }

    private static void genericsTypeNames(GenericsType[] genericsTypes, StringBuilder builder) {
        for( int i = 0; i < genericsTypes.length; i++ ) {
            if( i > 0 )
                builder.append(", ");
            builder.append(getName(genericsTypes[i].getType()));
        }
    }

    public static String getName(Class type) {
        return getName(type.getSimpleName());
    }

    public static String getName(String name) {
        if( "Object".equals(name) )
            return "?";
        return name;
    }

}
