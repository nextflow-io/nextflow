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
package nextflow.config.spec;

import java.lang.reflect.AnnotatedElement;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.lang.reflect.ParameterizedType;
import java.util.ArrayList;
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
 * Models a config spec.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public sealed interface SpecNode {
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

    private static SpecNode nextflowScope() {
        var enableOpts = new HashMap<String, SpecNode>();
        var previewOpts = new HashMap<String, SpecNode>();
        for( var field : FeatureFlagDsl.class.getDeclaredFields() ) {
            var fqName = field.getAnnotation(FeatureFlag.class).value();
            var names = fqName.split("\\.");
            var simpleName = names[names.length - 1];
            var desc = annotatedDescription(field, "");
            if( fqName.startsWith("nextflow.enable.") )
                enableOpts.put(simpleName, new Option(desc, optionTypes(field)));
            else if( fqName.startsWith("nextflow.preview.") )
                previewOpts.put(simpleName, new Option(desc, optionTypes(field)));
            else
                throw new IllegalArgumentException();
        }
        return new Scope(
            "",
            Map.ofEntries(
                Map.entry("enable", (SpecNode) new Scope("", enableOpts)),
                Map.entry("preview", (SpecNode) new Scope("", previewOpts))
            )
        );
    }

    /**
     * Initialize the `process` config scope from the set of
     * process directives.
     *
     * Directives with multiple method overloads are treated as
     * options with multiple supported types. Method overloads with
     * multiple parameters are ignored because they are not supported
     * in the configuration.
     */
    private static SpecNode processScope() {
        var description = """
            The `process` scope allows you to specify default directives for processes in your pipeline.
        
            [Read more](https://nextflow.io/docs/latest/config.html#process-configuration)
        """;
        var children = new HashMap<String, SpecNode>();
        for( var method : ProcessDsl.DirectiveDsl.class.getDeclaredMethods() ) {
            if( method.getParameters().length != 1 )
                continue;
            if( !children.containsKey(method.getName()) ) {
                var desc = annotatedDescription(method, "");
                children.put(method.getName(), new Option(desc, new ArrayList<>()));
            }
            var option = (Option) children.get(method.getName());
            var paramType = method.getParameterTypes()[0];
            option.types.add(paramType);
        }
        return new Scope(description, children);
    }
    
    private static String annotatedDescription(AnnotatedElement el, String defaultValue) {
        var annot = el.getAnnotation(Description.class);
        return annot != null ? annot.value() : defaultValue;
    }

    private static List<Class> optionTypes(Field field) {
        var result = new ArrayList<Class>();
        // use the field type
        result.add(field.getType());
        // append types from ConfigOption annotation if specified
        var annot = field.getAnnotation(ConfigOption.class);
        if( annot != null ) {
            for( var type : annot.types() )
                result.add(type);
        }
        return result;
    }

    /**
     * Models a config option that is defined through a DSL
     * instead of an assignment (i.e. `plugins`).
     */
    public static record DslOption(
        String description,
        Class dsl
    ) implements SpecNode {}

    /**
     * Models a config option.
     */
    public static record Option(
        String description,
        List<Class> types
    ) implements SpecNode {}

    /**
     * Models a config scope that contains custom named scopes
     * (e.g. `azure.batch.pools.<name>`).
     */
    public static record Placeholder(
        String description,
        String placeholderName,
        Scope scope
    ) implements SpecNode {}

    /**
     * Models a config scope.
     */
    public static record Scope(
        String description,
        Map<String, SpecNode> children
    ) implements SpecNode {
    
        /**
         * Get the spec node at the given path.
         *
         * @param names
         */
        public SpecNode getChild(List<String> names) {
            SpecNode node = this;
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
            var children = new HashMap<String, SpecNode>();
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
                        children.put(name, new Option(desc, optionTypes(field)));
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
