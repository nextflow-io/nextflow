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
package nextflow.script.types;

import java.nio.file.Path;
import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import groovy.lang.GString;
import groovy.lang.Tuple2;
import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.types.shim.ShimType;
import org.codehaus.groovy.GroovyBugError;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.GenericsType;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.tools.GenericsUtils;

/**
 * Utility constants and functions for working with Nextflow types.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class Types {

    public static final List<ClassNode> DEFAULT_IMPORTS = List.of(
        new ClassNode(Channel.class),
        new ClassNode(Duration.class),
        new ClassNode(MemoryUnit.class),
        new ClassNode(Path.class)
    );

    /**
     * Determine whether a method has a non-void return type.
     *
     * @param node
     */
    public static boolean hasReturnType(MethodNode node) {
        var returnType = node.getReturnType();
        if( returnType.isGenericsPlaceHolder() )
            return true;
        return !ClassHelper.isObjectType(returnType) && !ClassHelper.isPrimitiveVoid(returnType);
    }

    /**
     * Determine whether a source type can be assigned to a target type.
     *
     * @see StaticTypeCheckingSupport.checkCompatibleAssignmentTypes()
     *
     * @param target
     * @param source
     */
    public static boolean isAssignableFrom(ClassNode target, ClassNode source) {
        if( ClassHelper.isObjectType(target) )
            return true;
        if( target.equals(source) )
            return true;
        return isAssignableFrom(target.getTypeClass(), source.getTypeClass());
    }

    public static boolean isAssignableFrom(Class target, Class source) {
        target = normalize(target);
        source = normalize(source);
        if( target == Integer.class && source == Number.class )
            return false;
        if( target == Number.class && source == Integer.class )
            return true;
        return target.equals(source);
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

        if( type.isDerivedFrom(ClassHelper.TUPLE_TYPE) )
            return tupleName(type);

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
        else if( hasTypeClass(type) )
            builder.append(getName(type.getTypeClass()));
        else
            builder.append(getName(type.getNameWithoutPackage()));

        if( !placeholder && type.isUsingGenerics() ) {
            builder.append('<');
            genericsTypeNames(type.getGenericsTypes(), builder);
            builder.append('>');
        }

        return builder.toString();
    }

    private static void genericsTypeNames(GenericsType[] genericsTypes, StringBuilder builder) {
        for( int i = 0; i < genericsTypes.length; i++ ) {
            if( i > 0 )
                builder.append(", ");
            builder.append(getName(genericsTypes[i].getType()));
        }
    }

    private static boolean hasTypeClass(ClassNode type) {
        try {
            type.getTypeClass();
            return true;
        }
        catch( GroovyBugError e ) {
            return false;
        }
    }

    public static String getName(Class type) {
        return getName(normalize(type).getSimpleName());
    }

    public static String getName(String name) {
        if( "Object".equals(name) )
            return "?";
        return name;
    }

    private static final List<Class> STANDARD_TYPES = List.of(
        Boolean.class,
        Duration.class,
        Integer.class,
        List.class,
        Map.class,
        MemoryUnit.class,
        Number.class,
        Path.class,
        Set.class,
        String.class,
        Tuple2.class
    );

    private static final Map<Class,Class> PRIMITIVE_TYPES = Map.ofEntries(
        Map.entry(boolean.class, Boolean.class),
        Map.entry(double.class,  Number.class),
        Map.entry(float.class,   Number.class),
        Map.entry(int.class,     Integer.class),
        Map.entry(long.class,    Integer.class),
        Map.entry(GString.class, String.class)
    );

    /**
     * Determine the canonical Nextflow type for a given class.
     *
     * @param type
     */
    private static Class normalize(Class type) {
        if( type.isPrimitive() && PRIMITIVE_TYPES.containsKey(type) )
            return PRIMITIVE_TYPES.get(type);
        var queue = new LinkedList<Class>();
        queue.add(type);
        while( !queue.isEmpty() ) {
            var c = queue.remove();
            if( c == null )
                continue;
            if( STANDARD_TYPES.contains(c) )
                return c;
            queue.add(c.getSuperclass());
            for( var ic : c.getInterfaces() )
                queue.add(ic);
        }
        return type;
    }

    /**
     * Mapping of Java standard types to "shim" types which
     * provide method signatures and documentation.
     */
    public static final Map<ClassNode,ClassNode> SHIM_TYPES = shimTypes();

    private static Map<ClassNode,ClassNode> shimTypes() {
        var types = List.of(
            nextflow.script.types.shim.Bag.class,
            nextflow.script.types.shim.List.class,
            nextflow.script.types.shim.Map.class,
            nextflow.script.types.shim.Path.class,
            nextflow.script.types.shim.Set.class,
            nextflow.script.types.shim.String.class
        );
        return types.stream().collect(Collectors.toMap(
            (clazz) -> {
                var shim = clazz.getAnnotation(ShimType.class).value();
                return ClassHelper.makeCached(shim);
            },
            (clazz) -> ClassHelper.makeCached(clazz)
        ));
    }

}
