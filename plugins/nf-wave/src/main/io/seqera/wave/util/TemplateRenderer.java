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
 *
 */

package io.seqera.wave.util;

import java.io.InputStream;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Template rendering helper
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class TemplateRenderer {

    private static final Pattern PATTERN = Pattern.compile("\\{\\{([^}]+)}}");

    private static final Pattern VAR1 = Pattern.compile("(\\s*)\\{\\{([\\d\\w_-]+)}}(\\s*$)");
    private static final Pattern VAR2 = Pattern.compile("(?<!\\$)\\{\\{([\\d\\w_-]+)}}");

    private List<String> ignoreNames = List.<String>of();

    public TemplateRenderer withIgnore(String... names) {
        return withIgnore(List.of(names));
    }

    public TemplateRenderer withIgnore(List<String> names) {
        if( names!=null ) {
            ignoreNames = List.copyOf(names);
        }
        return this;
    }

    public String render(InputStream template, Map<String, String> binding) {
        String str = new Scanner(template).useDelimiter("\\A").next();
        return render(str, binding);
    }

    public String render(String template, Map<String, String> binding) {
        final String[] lines = template.split("(?<=\n)");
        final StringBuilder result = new StringBuilder();
        for( String it : lines ) {
            if( it==null || it.startsWith("##"))
                continue;
            final String resolved = replace0(it, binding);
            if( resolved!=null )
                result.append(resolved);
        }
        return result.toString();
    }

    /**
     * Simple template helper class replacing all variable enclosed by {{..}}
     * with the corresponding variable specified in a map object
     *
     * @param template The template string
     * @param binding The binding {@link Map}
     * @return The templated having the variables replaced with the corresponding value
     */
     String replace1(CharSequence template, Map<String, ?> binding) {
        Matcher matcher = PATTERN.matcher(template);

        // Loop through each matched variable placeholder
        StringBuilder builder = new StringBuilder();
        boolean isNull=false;
        while (matcher.find()) {
            String variable = matcher.group(1);

            // Check if the variable exists in the values map
            if (binding.containsKey(variable)) {
                Object value = binding.get(variable);
                String str = value!=null ? value.toString() : "";
                isNull |= value==null;
                matcher.appendReplacement(builder, str);
            }
            else if( !ignoreNames.contains(variable) ) {
                throw new IllegalArgumentException(String.format("Unable to resolve template variable: {{%s}}", variable));
            }
        }
        matcher.appendTail(builder);

        final String result = builder.toString();
        return !isNull || !result.isBlank() ? result : null;
    }

    String replace0(String line, Map<String,String> binding) {
        if( line==null || line.length()==0 )
            return line;

        Matcher matcher = VAR1.matcher(line);
        if( matcher.matches() ) {
            final String name = matcher.group(2);
            if( ignoreNames.contains(name) )
                return line;
            if( !binding.containsKey(name) )
                throw new IllegalArgumentException("Missing template key: "+name);
            final String prefix = matcher.group(1);
            final String value = binding.get(name);
            if( value==null )
                return null; // <-- return null to skip this line

            final StringBuilder result = new StringBuilder();
            final String[] multi = value.split("(?<=\n)");
            for (String s : multi) {
                result.append(prefix);
                result.append(s);
            }
            result.append( matcher.group(3) );
            return result.toString();
        }

        final StringBuilder result = new StringBuilder();
        while( (matcher=VAR2.matcher(line)).find() ) {
            final String name = matcher.group(1);
            if( !binding.containsKey(name) && !ignoreNames.contains(name)) {
                throw new IllegalArgumentException("Missing template key: "+name);
            }
            final String value = !ignoreNames.contains(name)
                    ? (binding.get(name)!=null ? binding.get(name) : "")
                    : "{{"+name+"}}";
            final int p = matcher.start(1);
            final int q = matcher.end(1);

            result.append(line.substring(0,p-2));
            result.append(value);
            line = line.substring(q+2);
        }
        result.append(line);
        return result.toString();
    }

}
