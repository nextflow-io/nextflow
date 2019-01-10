/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.util

import java.nio.file.Files
import java.nio.file.Path
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
/**
 * A basic mustache like parser
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MustacheTemplateEngine {

    @PackageScope static Pattern VAR1 = ~/(\s*)\{\{([\d\w_-]+)}}(\s*)$/

    @PackageScope static Pattern VAR2 = ~/(?<!\$)\{\{([\d\w_-]+)}}/

    String render(Path template, Map<String,String> binding) {
        def reader = Files.newBufferedReader(template)
        try {
            render(reader,binding)
        }
        finally {
            reader.close()
        }
    }

    protected boolean accept(String line) {
        return true
    }

    String render(String template, Map<String,String> binding) {
        render(new StringReader(template), binding)
    }

    String render(Reader template, Map<String,String> binding) {
        final result = new StringBuilder()

        final reader = template instanceof BufferedReader ? (BufferedReader)template : new BufferedReader(template)
        String line
        while( (line=reader.readLine()) != null ) {
            if( !accept(line) )
                continue
            def newLine = replace0(line, binding)
            if( newLine == null )
                continue
            result.append(newLine)
            if( !newLine.endsWith('\n') )
                result.append('\n')
        }

        return result.toString()
    }

    @PackageScope
    String replace0(String line, Map<String,String> binding) {
        if( !line )
            return line

        def matcher = VAR1.matcher(line)
        if( matcher.matches() ) {
            final name = matcher.group(2)
            if( !binding.containsKey(name) )
                throw new IllegalArgumentException("Missing template key: $name")
            final String prefix = matcher.group(1)
            final String value = binding.get(name)
            if( !value )
                return null // <-- return null to skip this line

            def result = new StringBuilder()
            final multi = value.readLines()
            for( int i=0; i<multi.size(); i++ ) {
                if(i) result.append('\n')
                result.append(prefix)
                result.append(multi[i])
            }
            return result
        }

        def result = new StringBuilder()
        while( (matcher=VAR2.matcher(line)).find() ) {
            def name = matcher.group(1)
            if( !binding.containsKey(name))
                throw new IllegalArgumentException("Missing template key: $name")
            def value = binding.get(name)
            def p = matcher.start(1)
            def q = matcher.end(1)

            result.append(line.substring(0,p-2))
            result.append(value)
            line = line.substring(q+2)
        }
        result.append(line)
        return result
    }

}
