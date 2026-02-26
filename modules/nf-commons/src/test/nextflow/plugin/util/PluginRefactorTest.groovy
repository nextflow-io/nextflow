/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.plugin.util


import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PluginRefactorTest extends Specification {


    def "should normalize strings into PascalCase class names"() {
        given:
        def refactor = new PluginRefactor()
        expect:
        refactor.normalizeToClassName(INPUT) == EXPECTED

        where:
        INPUT                           || EXPECTED
        "  my-cool_class name!  "       || "MyCoolClassName"
        "acme's #1 utility!"            || "AcmeS1Utility"
        "hello_world"                   || "HelloWorld"
        "123 start here"                || "123StartHere"
        "   mixed---separators___"      || "MixedSeparators"
        "alreadyPascalCase"             || "AlreadyPascalCase"
        "foo-plugin"                    || "Foo"
        "nf-hello"                      || "Hello"
        "nf-my-plugin"                  || "My"
        "NfHelloWorld"                  || "HelloWorld"
        "my-nf-plugin"                  || "MyNf"
    }

    def "should normalize strings into kebab-case names"() {
        given:
        def refactor = new PluginRefactor()
        expect:
        refactor.normalizeToKebabCase(INPUT) == EXPECTED

        where:
        INPUT                           || EXPECTED
        "MyCoolClassName"               || "my-cool-class-name"
        "my-cool_class name!"           || "my-cool-class-name"
        "acme's #1 utility!"            || "acme-s-1-utility"
        "hello_world"                   || "hello-world"
        "123 start here"                || "123-start-here"
        "Mixed---separators___"         || "mixed-separators"
        "AlreadyKebab-Case"             || "already-kebab-case"
        "CamelCaseInput"                || "camel-case-input"
        "HTMLParserUtility"             || "html-parser-utility"
    }


    def "should normalize strings into valid Java package name segments"() {
        given:
        def refactor = new PluginRefactor()
        expect:
        refactor.normalizeToPackageNameSegment(INPUT) == EXPECTED

        where:
        INPUT                   || EXPECTED
        "Acme Tools!"           || "acmetools"
        "123-Utility"           || "utility"
        "com.example"           || "comexample"
        "__System_Core__"       || "systemcore"
        "1st-Package"           || "stpackage"
        "My.App.Service"        || "myappservice"
        "!@#\$%^&*()"           || null
        ""                      || null
    }
}
