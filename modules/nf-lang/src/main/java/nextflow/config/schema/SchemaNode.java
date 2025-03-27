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
package nextflow.config.schema;

import java.lang.reflect.AnnotatedElement;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.lang.reflect.ParameterizedType;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import nextflow.config.scopes.Config;
import nextflow.script.dsl.Description;
import nextflow.script.dsl.DslScope;
import nextflow.script.dsl.FeatureFlag;
import nextflow.script.dsl.FeatureFlagDsl;
import nextflow.script.dsl.ProcessDsl;

/**
 * Models the config schema.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public sealed interface SchemaNode {
    String description();

    public static final Scope ROOT = rootScope();

    private static Scope rootScope() {
        var result = Scope.of(Config.class, "");
        // derive `nextflow` config options from feature flags.
        result.children().put("nextflow", nextflowScope());
        // derive `process` config options from process directives.
        result.children().put("process", processScope());
        return result;
    }

    private static SchemaNode nextflowScope() {
        var enableOpts = new HashMap<String, SchemaNode>();
        var previewOpts = new HashMap<String, SchemaNode>();
        for( var field : FeatureFlagDsl.class.getDeclaredFields() ) {
            var fqName = field.getAnnotation(FeatureFlag.class).value();
            var names = fqName.split("\\.");
            var simpleName = names[names.length - 1];
            var desc = annotatedDescription(field, "");
            if( fqName.startsWith("nextflow.enable.") )
                enableOpts.put(simpleName, new Option(desc, optionType(field)));
            else if( fqName.startsWith("nextflow.preview.") )
                previewOpts.put(simpleName, new Option(desc, optionType(field)));
            else
                throw new IllegalArgumentException();
        }
        return new Scope(
            "",
            Map.ofEntries(
                Map.entry("enable", (SchemaNode) new Scope("", enableOpts)),
                Map.entry("preview", (SchemaNode) new Scope("", previewOpts))
            )
        );
    }

    private static SchemaNode processScope() {
        var description = """
            The `process` scope allows you to specify default directives for processes in your pipeline.
        
            [Read more](https://nextflow.io/docs/latest/config.html#process-configuration)
        """;
        var children = new HashMap<String, SchemaNode>();
        for( var method : ProcessDsl.DirectiveDsl.class.getDeclaredMethods() ) {
            var desc = annotatedDescription(method, "");
            children.put(method.getName(), new Option(desc, optionType(method)));
        }
        return new Scope(description, children);
    }
    
    private static String annotatedDescription(AnnotatedElement el, String defaultValue) {
        var annot = el.getAnnotation(Description.class);
        return annot != null ? annot.value() : defaultValue;
    }

    private static Class optionType(AnnotatedElement element) {
        if( element instanceof Field field ) {
            return field.getType();
        }
        if( element instanceof Method method ) {
            // use the return type if config option is not a directive
            var returnType = method.getReturnType();
            if( returnType != void.class )
                return returnType;
            // other use the type of the last parameter
            var paramTypes = method.getParameterTypes();
            if( paramTypes.length > 0 )
                return paramTypes[paramTypes.length - 1];
        }
        return null;
    }

    /**
     * Models a config option that is defined through a DSL
     * instead of an assignment (i.e. `plugins`).
     */
    public static record DslOption(
        String description,
        Class dsl
    ) implements SchemaNode {}

    /**
     * Models a config option.
     */
    public static record Option(
        String description,
        Class type
    ) implements SchemaNode {}

    /**
     * Models a config scope that contains custom named scopes
     * (e.g. `azure.batch.pools.<name>`).
     */
    public static record Placeholder(
        String description,
        String placeholderName,
        Scope scope
    ) implements SchemaNode {}

    /**
     * Models a config scope.
     */
    public static record Scope(
        String description,
        Map<String, SchemaNode> children
    ) implements SchemaNode {
    
        /**
         * Get the schema node at the given path.
         *
         * @param names
         */
        public SchemaNode getChild(List<String> names) {
            SchemaNode node = this;
            for( var name : names ) {
                if( node instanceof Scope sn )
                    node = sn.children().get(name);
                else if( node instanceof Placeholder pn )
                    node = pn.scope();
                else
                    return null;
            }
            return node;
        }
    
        /**
         * Get the config dsl option at the given path.
         *
         * @param names
         */
        public DslOption getDslOption(List<String> names) {
            return getChild(names) instanceof DslOption option ? option : null;
        }
    
        /**
         * Get the config option at the given path.
         *
         * @param names
         */
        public Option getOption(List<String> names) {
            return getChild(names) instanceof Option option ? option : null;
        }

        /**
         * Get the config scope at the given path.
         *
         * @param names
         */
        public Scope getScope(List<String> names) {
            return getChild(names) instanceof Scope scope ? scope : null;
        }
    
        /**
         * Create a scope node from a ConfigScope class.
         *
         * @param scope
         * @param description
         */
        public static Scope of(Class<? extends ConfigScope> scope, String description) {
            var children = new HashMap<String, SchemaNode>();
            for( var field : scope.getDeclaredFields() ) {
                var name = field.getName();
                var type = field.getType();
                var desc = annotatedDescription(field, description);
                var placeholderName = field.getAnnotation(PlaceholderName.class);
                // fields annotated with @ConfigOption are config options
                if( field.getAnnotation(ConfigOption.class) != null ) {
                    if( DslScope.class.isAssignableFrom(type) )
                        children.put(name, new DslOption(desc, type));
                    else
                        children.put(name, new Option(desc, optionType(field)));
                }
                // fields of type ConfigScope are nested config scopes
                else if( ConfigScope.class.isAssignableFrom(type) ) {
                    children.put(name, Scope.of((Class<? extends ConfigScope>) type, desc));
                }
                // fields of type Map<String, ConfigScope> are placeholder scopes
                else if( Map.class.isAssignableFrom(type) && placeholderName != null ) {
                    var pt = (ParameterizedType)field.getGenericType();
                    var valueType = (Class<? extends ConfigScope>)pt.getActualTypeArguments()[1];
                    children.put(name, new Placeholder(desc, placeholderName.value(), Scope.of(valueType, desc)));
                }
            }
            return new Scope(description, children);
        }
    }

}
