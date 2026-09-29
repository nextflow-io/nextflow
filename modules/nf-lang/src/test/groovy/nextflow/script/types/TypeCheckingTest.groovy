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

package nextflow.script.types

import nextflow.script.ast.ASTNodeMarker
import nextflow.script.ast.FeatureFlagNode
import nextflow.script.control.ScriptParser
import nextflow.script.control.ScriptResolveVisitor
import nextflow.script.control.TypeCheckingVisitor
import nextflow.script.dsl.Types
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassHelper
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.messages.SyntaxErrorMessage
import org.codehaus.groovy.syntax.SyntaxException
import spock.lang.Ignore
import spock.lang.Specification
import spock.lang.Unroll

import static nextflow.script.types.TypeCheckingUtils.*

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypeCheckingTest extends Specification {

    SourceUnit parse(String contents) {
        def scriptParser = new ScriptParser()
        def source = scriptParser.parse('main.nf', contents.stripIndent())
        assert source.getAST() != null
        source.getAST().addFeatureFlag(new FeatureFlagNode("nextflow.enable.types", new ConstantExpression(true)))
        new ScriptResolveVisitor(source, scriptParser.compiler().compilationUnit(), Types.DEFAULT_SCRIPT_IMPORTS, Collections.emptyList()).visit()
        new TypeCheckingVisitor(source).visit()
        return source
    }

    Statement parseStatement(String contents) {
        def source = parse(contents)
        assert !source.getErrorCollector().hasErrors()
        def entry = source.getAST().getEntry()
        assert entry.main instanceof BlockStatement
        return entry.main.statements.last()
    }

    Expression parseExpression(String contents) {
        def stmt = parseStatement(contents)
        assert stmt instanceof ExpressionStatement
        return stmt.expression
    }

    def getErrors(String contents) {
        final source = parse(contents)
        final errorCollector = source.getErrorCollector()
        if( !errorCollector.hasErrors() )
            return Collections.emptyList()
        return errorCollector.getErrors().stream()
            .filter(e -> e instanceof SyntaxErrorMessage)
            .map(e -> e.cause)
            .sorted(ERROR_COMPARATOR)
            .toList()
    }

    static final Comparator<SyntaxException> ERROR_COMPARATOR = (SyntaxException a, SyntaxException b) -> {
        return a.getStartLine() != b.getStartLine()
            ? a.getStartLine() - b.getStartLine()
            : a.getStartColumn() - b.getStartColumn()
    }

    boolean check(String contents, String message) {
        final errors = getErrors(contents)
        if( message != null ) {
            assert errors.size() == 1
            assert errors[0].getOriginalMessage() == message
        }
        else {
            assert errors.size() == 0
        }
        return true
    }

    boolean checkType(ASTNode exp, Class type) {
        assert Types.isEqual(getType(exp), ClassHelper.makeCached(type))
        return true
    }

    def 'should warn about a process `when` section' () {
        when:
        def errors = getErrors(
            '''\
            nextflow.enable.types = true

            process hello {
                when:
                task.ext.when

                exec:
                println 'hello!'
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 4
        errors[0].getStartColumn() == 5
        errors[0].isSoftError()
        errors[0].getOriginalMessage() == "Process `when` section is discouraged with static typing -- use conditional logic in the calling workflow instead"
    }

    @Unroll
    def 'should report legacy type annotations' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE                                          | ERROR
        "String x = 'a'"                                | "Legacy type annotations are not supported with static typing -- use `name: Type` instead"
        "def String x = 'a'"                            | "Legacy type annotations are not supported with static typing -- use `name: Type` instead"
        "def hello(String x) { return x }"              | "Legacy type annotations are not supported with static typing -- use `name: Type` instead"
        "def x: String = 'a'"                           | null
        "def hello(x: String) { return x }"             | null
    }

    @Unroll
    def 'should resolve a static method call' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE                                          | ERROR
        "String.format('%d', 1)"                        | null
        "Integer.parseInt('1')"                         | null
        "Math.max(1, 2)"                                | null
        "groovy.json.JsonOutput.toJson([1, 2])"         | null
        "String.nope()"                                 | "Unrecognized method `nope` for type String"
    }

    def 'should check a parameter declaration' () {
        when:
        def errors = getErrors(
            '''\
            params {
                input: String = 3.14
            }

            workflow {
                params
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 2
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "Parameter 'input' with type String cannot be assigned to default value with type Float"

        when:
        def exp = parseExpression(
            '''\
            params {
                input: String = 'input'
            }

            workflow {
                params
            }
            '''
        )
        def cn = getType(exp)
        then:
        checkType(exp, ParamsMap)
        cn.getField('input') != null
    }

    def 'should check a workflow emit' () {
        when:
        def errors = getErrors(
            '''\
            nextflow.enable.types = true

            workflow hello {
                emit:
                a: String = 42
                b: Integer = 1
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 5
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "Assignment target with type String cannot be assigned to value with type Integer"

        when:
        errors = getErrors(
            '''\
            nextflow.enable.types = true

            workflow hello {
                emit:
                a: Integer = 42
                b: Integer = 1
            }
            '''
        )
        then:
        errors.size() == 0
    }

    def 'should warn about a single named output' () {
        expect:
        check(
            '''\
            nextflow.enable.types = true

            workflow hello {
                emit:
                result: Integer = 42
            }
            ''',
            "Name should be omitted for a single emit"
        )
        check(
            '''\
            nextflow.enable.types = true

            process hello {
                output:
                sample: String = 'hi'

                script:
                ''
            }
            ''',
            "Name should be omitted for a single output"
        )
    }

    def 'should check output path directive' () {
        when:
        def errors = getErrors(
            '''\
            nextflow.enable.types = true

            workflow {
                main:
                ch_samples = channel.empty()

                publish:
                samples = ch_samples
            }

            output {
                samples: Channel<Sample> {
                    path { s ->
                        s.id >> 42
                        s.fastq >> params.save_fastq ? 'fastq' : null
                    }
                }
            }

            record Sample {
                id: String
                fastq: Path
            }
            '''
        )
        then:
        errors.size() == 3
        errors[0].getStartLine() == 14
        errors[0].getStartColumn() == 13
        errors[0].getOriginalMessage() == "Publish source should be a Path or Iterable<Path> but was specified as a String"
        errors[1].getStartLine() == 14
        errors[1].getStartColumn() == 21
        errors[1].getOriginalMessage() == "Publish target should be a String but was specified as a Integer"
        errors[2].getStartLine() == 15
        errors[2].getStartColumn() == 13
        errors[2].getOriginalMessage().contains "Statement is not a valid publish statement"

        when:
        errors = getErrors(
            '''\
            nextflow.enable.types = true

            workflow {
                main:
                ch_samples = channel.empty()

                publish:
                samples = ch_samples
            }

            output {
                samples: Channel<Sample> {
                    path { s ->
                        s.fastq >> (params.save_fastq ? 'fastq' : null)
                    }
                }
            }

            record Sample {
                id: String
                fastq: Path
            }
            '''
        )
        then:
        errors.size() == 0
    }

    def 'should infer the return type of an untyped function' () {
        when:
        def source = parse(
            '''\
            def make() {
                return channel.of(1, 2, 3)
            }

            workflow {
                make()
            }
            '''
        )
        def fn = source.getAST().getFunctions().first()
        def call = source.getAST().getEntry().main.statements.last().expression
        then:
        Types.getName(getType(call)) == 'Channel<Integer>'
        // the function itself stays untyped, since code generation reads it
        ClassHelper.isDynamicTyped(fn.getReturnType())
    }

    def 'should check a return statement' () {
        when:
        def errors = getErrors(
            '''\
            def hello() -> String {
                return 42
            }
            '''
        )
        then:
        errors.size() == 1
        errors[0].getStartLine() == 2
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "Return value with type Integer does not match the declared return type (String)"

        when:
        errors = getErrors(
            '''\
            def hello() -> String {
                return 'Hello!'
            }
            '''
        )
        then:
        errors.size() == 0

        when: 'the return statement is implicit'
        errors = getErrors(
            '''\
            def hello() -> String {
                println 'Hello!'
            }
            '''
        )
        then: 'the error is reported against the statement it was inferred from'
        errors.size() == 1
        errors[0].getStartLine() == 2
        errors[0].getStartColumn() == 5
        errors[0].getOriginalMessage() == "Return value with type void does not match the declared return type (String)"
    }

    def 'should check an assignment' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE                              | ERROR
        'def x: String ; x = 42'            | "Assignment target with type String cannot be assigned to value with type Integer"
        "def x: List<String> ; x = [42]"    | "Assignment target with type List<String> cannot be assigned to value with type List<Integer>"
        "def x: List<String> ; x = []"      | null
        "def ch: Channel<Path> = channel.empty()" | null
    }

    @Unroll
    def 'should check a compound assignment' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE                          | ERROR
        'def x = 1 ; x += 1'            | null
        'def x = 1 ; x -= 1'            | null
        'def x = 1 ; x *= 2'            | null
        'def x = 1 ; x <<= 2'           | null
        "def s = 'a' ; s += 'b'"        | null
        'def d = 1.h ; d *= 2.0'        | null
        "def s = 'a' ; s *= 'b'"        | "The `*=` operator is not defined for operands with types String and String"
        "def x = 1 ; x += 'a'"          | "The `+=` operator is not defined for operands with types Integer and String"
    }

    def 'should infer the type of a variable declaration or assignment' () {
        when:
        def exp = parseExpression(
            '''
            def x = 'hello'
            '''
        )
        def target = exp.getLeftExpression()
        then:
        checkType(target, String)

        when:
        exp = parseExpression(
            '''
            def x
            x = 'hello'
            '''
        )
        target = exp.getLeftExpression()
        then:
        checkType(target, String)

        when:
        exp = parseExpression(
            '''
            x = null
            x = 'hello'
            '''
        )
        target = exp.getLeftExpression()
        then:
        checkType(target, String)

        when:
        exp = parseExpression(
            '''
            def x = 1
            x += 1
            '''
        )
        target = exp.getLeftExpression()
        then:
        checkType(target, Integer)
    }

    @Unroll
    def 'should check a function call' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE          | ERROR
        'env()'         | "Function `env` expects 1 argument(s) but received 0"
        'env(42)'       | "Argument with type Integer is not compatible with parameter of type String"
        "env('HELLO')"  | null
    }

    @Unroll
    def 'should check a function call with overloads' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE              | ERROR
        'file()'            | "Function `file` (with multiple signatures) was called with incorrect number of arguments and/or incorrect argument types"
        'file(42)'          | "Function `file` (with multiple signatures) was called with incorrect number of arguments and/or incorrect argument types"
        "file('input.txt')" | null
        "file('input.txt', checkIfExists: true)" | null
    }

    @Unroll
    def 'should check a function call with named params' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE                                      | ERROR
        "file('input.txt', checkIfExist: true)"     | "Named param `checkIfExist` is not defined"
        "file('input.txt', checkIfExists: 'true')"  | "Named param `checkIfExists` expects a Boolean but received a String"
        "file('input.txt', checkIfExists: true)"    | null
    }

    def 'should check a function call with record parameters' () {
        expect:
        check(
            '''\
            def hello(s: Sample) {
            }

            record Sample {
                id: String
                fastq: Path
            }

            workflow {
                hello( record(id: '1') )
            }
            ''',
            'Argument with type Record {\n    id: String\n} is not compatible with parameter of type Sample'
        )
        and:
        check(
            '''\
            def hello(s: Sample) {
            }

            record Sample {
                id: String
                fastq: Path
            }

            workflow {
                hello( record(id: '1', fastq: file('1.fastq')) )
            }
            ''',
            null
        )
    }

    @Unroll
    def 'should check call arguments that infer conflicting type arguments' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE                          | ERROR
        "channel.of(1, 'a')"            | "Argument with type String is not compatible with parameter of type Integer"
        "channel.of('a', 1)"            | "Argument with type Integer is not compatible with parameter of type String"
        "channel.of(1, 2, 'a')"         | "Argument with type String is not compatible with parameter of type Integer"
        "channel.of([1], ['a'])"        | "Argument with type List<String> is not compatible with parameter of type List<Integer>"
        "channel.of(1, 2)"              | null
        "channel.of([1], [2])"          | null
        "channel.of(1).reduce(0) { a, b -> 'x' }" | "Return value with type String does not match the declared return type (Integer)"
    }

    @Unroll
    def 'should check a namespaced function call' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE                  | ERROR
        'channel.queue()'       | "Unrecognized method `queue` for namespace `channel`"
        'channel.value()'       | "Function `value` expects 1 argument(s) but received 0"
        'channel.of()'          | null
        'channel.of(1, 2, 3)'   | null
    }

    @Unroll
    def 'should check a method call' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE                                  | ERROR
        'workflow.outputDir.name()'             | "Unrecognized method `name` for type Path"
        'workflow.outputDir.resolve()'          | "Function `resolve` expects 1 argument(s) but received 0"
        "workflow.outputDir.resolve('hello')"   | null
    }

    @Unroll
    def 'should check a binary expression' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE          | ERROR
        "2 + '2'"       | "The `+` operator is not defined for operands with types Integer and String"
        "2 + 2"         | null
        "2 == '2'"      | "The `==` operator is not defined for operands with types Integer and String"
        "2 == 2"        | null
        "2 in ['2']"    | "The `in` operator is not defined for operands with types Integer and List<String>"
        "2 in [2]"      | null
        "[2] == ['2']"  | "The `==` operator is not defined for operands with types List<Integer> and List<String>"
        "[2] == [2]"    | null
        "['a', 'b', 'c'][0]"        | null
        "['a', 'b', 'c'][0..-2]"    | null
        "['a', 'b', 'c'][0..1]"     | null
        "['a', 'b', 'c'][1..<3]"    | null
    }

    @Unroll
    def 'should resolve a binary expression' () {
        given:
        def exp = parseExpression(SOURCE)
        expect:
        Types.getName(getType(exp)) == TYPE

        where:
        SOURCE          | TYPE
        '1 <=> 2'       | 'Integer'
        "'a' <=> 'b'"   | 'Integer'
        '1 < 2'         | 'Boolean'
    }

    @Unroll
    def 'should resolve list slicing to the list type' () {
        given:
        def exp = parseExpression(SOURCE)
        expect:
        Types.getName(getType(exp)) == TYPE

        where:
        SOURCE                      | TYPE
        "['a', 'b', 'c'][0]"        | 'String'
        "['a', 'b', 'c'][0..-2]"    | 'List<String>'
    }

    def 'should recognize duration literals' () {
        expect:
        check(
            '''\
            (1.d + 4.h + 15.m + 30.s + 500.ms).toMillis()
            ''',
            null
        )
    }

    def 'should recognize memory unit literals' () {
        expect:
        check(
            '''\
            (1.GB + 500.MB + 4.KB + 64.B).toBytes()
            ''',
            null
        )
    }

    @Unroll
    def 'should check a ternary expression' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE                  | ERROR
        "true ? 42 : '42'"      | "Conditional expression has inconsistent types -- true branch has type Integer but false branch has type String"
        "true ? ['v'] : [:]"    | "Conditional expression has inconsistent types -- true branch has type List<String> but false branch has type Map"
        "true ? 42 : null"      | null
    }

    def 'should resolve a ternary expression' () {
        when:
        def exp = parseExpression(SOURCE)
        then:
        Types.getName(getType(exp)) == TYPE

        where:
        SOURCE                  | TYPE
        "true ? 42 : null"      | 'Integer?'
        "true ? [k: 'v'] : [:]" | 'Map<String, String>'
        "true ? [:] : [k: 'v']" | 'Map<String, String>'
        "true ? ['v'] : []"     | 'List<String>'
        "true ? [] : ['v']"     | 'List<String>'
    }

    def 'should check an elvis expression' () {
        when:
        def exp = parseExpression(
            '''
            workflow {
                def meta: Map<String,String> = [:]
                meta.id ?: 'id'
            }
            '''
        )
        then:
        checkType(exp, String)
    }

    @Unroll
    def 'should check a list expression' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE          | ERROR
        "[1, 2, '3']"   | "List expression has inconsistent element types -- some elements have type Integer while others have type String"
        "[1, 2, 3]"     | null
    }

    @Unroll
    def 'should resolve generic types' () {
        given:
        def exp = parseExpression(SOURCE)
        expect:
        Types.getName(getType(exp)) == TYPE

        where:
        SOURCE                              | TYPE
        '[1, 2, 3]'                         | 'List<Integer>'
        '[x: 1, y: 2]'                      | 'Map<String, Integer>'
        "[x: '1', y: 2]"                    | 'Map<String, ?>'
        "([x: '1'] + [y: '2']).keySet()"    | 'Set<String>'
        "([x: '1'] + [y: '2']).values()"    | 'Bag<String>'
        '[1, 2, 3].first()'                 | 'Integer'
        'channel.value("42")'               | 'Value<String>'
        'channel.of(1, 2, 3)'               | 'Channel<Integer>'
        '[1, 2, 3].collect { v -> v * 2 }'  | 'Iterable<Integer>'
        'channel.of(1).map { v -> v * 2 }'  | 'Channel<Integer>'
        'channel.fromList(1..10).flatMap { n -> 1..n }' | 'Channel<Integer>'
    }

    def 'should resolve generic functional parameter types' () {
        when:
        def exp = parseExpression(
            """\
            [1, 2, 3].collect { v -> v * 2 }
            """
        )
        def method = exp.getNodeMetaData(ASTNodeMarker.METHOD_TARGET)
        def param = method.getParameters().last()
        then:
        Types.getName(param.getType()) == '(Integer) -> Integer'

        when:
        exp = parseExpression(
            """\
            channel.of(1, 2, 3).map { v -> v * 2 }
            """
        )
        method = exp.getNodeMetaData(ASTNodeMarker.METHOD_TARGET)
        param = method.getParameters().last()
        then:
        Types.getName(param.getType()) == '(Integer) -> Integer'
    }

    def 'should resolve the closure variants of toSorted' () {
        when: 'a transform closure'
        def exp = parseExpression(
            """\
            [1, 2, 3].toSorted { v -> -v }
            """
        )
        def method = exp.getNodeMetaData(ASTNodeMarker.METHOD_TARGET)
        then:
        Types.getName(method.getParameters().last().getType()) == '(Integer) -> Integer'
        Types.getName(getType(exp)) == 'List<Integer>'

        when: 'a comparator closure'
        exp = parseExpression(
            """\
            [1, 2, 3].toSorted { a, b -> a <=> b }
            """
        )
        method = exp.getNodeMetaData(ASTNodeMarker.METHOD_TARGET)
        then:
        Types.getName(method.getParameters().last().getType()) == '(Integer, Integer) -> Integer'
        Types.getName(getType(exp)) == 'List<Integer>'
    }

    @Unroll
    def 'should check a unary expression' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE  | ERROR
        "-'42'" | "The `-` operator is not defined for an operand with type String"
        "-42"   | null
    }

    @Unroll
    def 'should check a cast expression' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE                  | ERROR
        '42 as List'            | "Value of type Integer cannot be cast to List"
        '[] as List<Path>'      | null
        "'24 h' as Duration"    | null
        "'8 GB' as MemoryUnit"  | null
    }

    def 'should check a record cast' () {
        expect:
        check(
            '''\
            workflow {
                record(id: '1', fastq_1: file('1_1.fastq'), fastq_2: file('1_2.fastq')) as FastqPair
            }

            record FastqPair {
                id: String
                fastq_1: Path
                fastq_2: Path
            }
            ''',
            null
        )
        check(
            '''\
            workflow {
                record(id: '1', fastq_1: file('1_1.fastq')) as FastqPair
            }

            record FastqPair {
                id: String
                fastq_1: Path
                fastq_2: Path
            }
            ''',
            'Record mismatch -- source record is missing field `fastq_2` required by FastqPair'
        )
        check(
            '''\
            workflow {
                record(id: '1', fastq_1: file('1_1.fastq'), fastq_2: '1_2.fastq') as FastqPair
            }

            record FastqPair {
                id: String
                fastq_1: Path
                fastq_2: Path
            }
            ''',
            'Record mismatch -- field `fastq_2` of FastqPair expects a Path but received a String'
        )
    }

    def 'should check a record channel assignment' () {
        expect:
        check(
            '''\
            workflow {
                def ch: Channel<Sample> = channel.of( record(id: '1', fastq: '1.fastq') )
            }

            record Sample {
                id: String
                fastq: String
            }
            ''',
            null
        )
        check(
            '''\
            workflow {
                def ch: Channel<Sample> = channel.of( record(id: '1') )
            }

            record Sample {
                id: String
                fastq: String
            }
            ''',
            'Assignment target with type Channel<Sample> cannot be assigned to value with type Channel<Record {\n    id: String\n}>'
        )
    }

    @Unroll
    def 'should check a property expression' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE                          | ERROR
        'workflow.output_dir'           | "Unrecognized property `output_dir` for namespace `workflow`"
        'workflow.outputDir'            | null
        'workflow.outputDir.fileName'   | "Unrecognized property `fileName` for type Path"
        'workflow.outputDir.name'       | null
    }

    def 'should check a map property' () {
        when:
        def exp = parseExpression(
            '''
            workflow {
                def meta: Map = [:]
                meta.id
            }
            '''
        )
        then:
        checkType(exp, Object)

        when:
        exp = parseExpression(
            '''
            workflow {
                def meta: Map = [:]
                meta = meta + [single_end: true]
                meta.id
            }
            '''
        )
        then:
        checkType(exp, Object)
    }

    def 'should check a workflow call' () {
        expect:
        check(
            '''
            nextflow.enable.types = true

            workflow hello {
                take:
                messages: Channel<String>

                main:
                messages.view()
            }

            workflow {
                hello()
            }
            ''',
            'Workflow `hello` expects 1 argument(s) but received 0'
        )
    }

    def 'should recognize workflow output type' () {
        when:
        def exp = parseExpression(
            '''\
            nextflow.enable.types = true

            workflow hello {
                emit:
                42
            }

            workflow {
                hello()
            }
            '''
        )
        def type = getType(exp)
        then:
        Types.getName(type) == 'Value<Integer>'

        when:
        exp = parseExpression(
            '''\
            nextflow.enable.types = true

            workflow hello {
                emit:
                foo: String = 'hello'
                bar: Integer = 42
            }

            workflow {
                hello().bar
            }
            '''
        )
        type = getType(exp)
        then:
        Types.getName(type) == 'Value<Integer>'
    }

    def 'should check a process call' () {
        expect:
        check(
            '''
            nextflow.enable.types = true

            process hello {
                input:
                message: String

                script:
                ''
            }

            workflow hello_flow {
                take:
                messages: Channel<String>

                main:
                hello( messages )
            }
            ''',
            null
        )
        and:
        check(
            '''
            nextflow.enable.types = true

            process hello {
                input:
                message: String

                script:
                ''
            }

            workflow {
                hello()
            }
            ''',
            'Process `hello` expects 1 argument(s) but received 0'
        )
        and:
        check(
            '''
            nextflow.enable.types = true

            process hello {
                input:
                message: String

                script:
                ''
            }

            workflow hello_flow {
                take:
                messages: Channel<Integer>

                main:
                hello( messages )
            }
            ''',
            'Argument with type Integer is not compatible with process input of type String'
        )
        and:
        check(
            '''
            nextflow.enable.types = true

            process hello {
                input:
                message: String
                target: String

                script:
                ''
            }

            workflow hello_flow {
                take:
                messages: Channel<String>
                targets: Channel<String>

                main:
                hello( messages, targets )
            }
            ''',
            'Process `hello` was called with multiple channel arguments which can lead to non-deterministic behavior -- make sure that at most one argument is a channel and that all other arguments are dataflow values'
        )
    }

    def 'should not check the return value of a closure with a void signature' () {
        expect: 'each() takes a Consumer, which discards the closure result'
        check(
            '''\
            def total = 0
            [1, 2, 3].each { n -> total = total + n }
            ''',
            null
        )
    }

    def 'should not check the return value of a closure with a boolean signature' () {
        expect: 'findAll() takes a Predicate, which coerces the closure result to a boolean'
        check(
            "['a', 'b', ''].findAll { s -> s }",
            null
        )
    }

    def 'should treat a wildcard type argument as unknown' () {
        expect: 'channel.topic() returns Channel<?>, so the element type is not known'
        check(
            "channel.topic('versions').view { v -> v.name }",
            null
        )
        and: 'the same call on a known element type is still checked'
        check(
            "channel.of('a').view { v -> v.name }",
            'Unrecognized property `name` for type String'
        )
    }

    @Unroll
    def 'should report a legacy operator as discouraged rather than unrecognized' () {
        expect:
        check("channel.of('a').${SOURCE}", "Operator `${OPERATOR}` is discouraged from use with static typing")

        where:
        SOURCE                          | OPERATOR
        "toSortedList { s -> s.size() }"| 'toSortedList'
        "toSortedList()"                | 'toSortedList'
        "count('a')"                    | 'count'
        "count()"                       | 'count'
        "max()"                         | 'max'
        "min()"                         | 'min'
        "sum()"                         | 'sum'
        "groupTuple()"                  | 'groupTuple'
        "collectFile()"                 | 'collectFile'
    }

    def 'should check an agent call' () {
        expect: 'an agent call is lifted and unwrapped like a process call'
        check(
            '''\
            nextflow.enable.types = true

            agent hello {
                input:
                message: String

                output:
                answer: String

                prompt:
                "say ${message}"
            }

            workflow {
                hello( channel.of('hi') ).view { v -> v.toUpperCase() }
            }
            ''',
            null
        )
        and:
        check(
            '''\
            nextflow.enable.types = true

            agent hello {
                input:
                message: String

                output:
                answer: String

                prompt:
                "say ${message}"
            }

            workflow {
                hello( channel.of(42) )
            }
            ''',
            'Argument with type Integer is not compatible with agent input of type String'
        )
    }

    def 'should check a process call with record inputs' () {
        expect:
        check(
            '''\
            nextflow.enable.types = true

            process hello {
                input:
                record(id: String, fastq: Path)

                exec:
                println '...'
            }

            workflow {
                hello( record(id: '1') )
            }
            ''',
            'Argument with type Record {\n    id: String\n} is not compatible with process input of type Record {\n    id: String\n    fastq: Path\n}'
        )
        and:
        check(
            '''\
            nextflow.enable.types = true

            process hello {
                input:
                record(id: String, fastq: Path)

                exec:
                println '...'
            }

            workflow {
                hello( record(id: '1', fastq: file('1.fastq')) )
            }
            ''',
            null
        )
    }

    def 'should recognize process output type' () {
        when:
        def exp = parseExpression(
            '''\
            nextflow.enable.types = true

            process hello {
                input:
                target: String

                output:
                "Hello, $target!"

                exec:
                true
            }

            workflow {
                hello( 'World' )
            }
            '''
        )
        def type = getType(exp)
        then:
        Types.getName(type) == 'Value<String>'

        when:
        exp = parseExpression(
            '''\
            nextflow.enable.types = true

            process hello {
                input:
                target: String

                output:
                tuple(target, "Hello, $target!")

                exec:
                true
            }

            workflow {
                hello( 'World' )
            }
            '''
        )
        type = getType(exp)
        then:
        Types.getName(type) == 'Value<Tuple<String, String>>'

        when:
        exp = parseExpression(
            '''\
            nextflow.enable.types = true

            process hello {
                input:
                target: String

                output:
                record(target: target, message: "Hello, $target!")

                exec:
                true
            }

            workflow {
                hello( 'World' )
            }
            '''
        )
        type = getType(exp)
        then:
        Types.getName(type) == 'Value<Record {\n    target: String\n    message: String\n}>'
    }

    def 'should report error for process .out property' () {
        expect:
        check(
            '''
            nextflow.enable.types = true

            process hello {
                output:
                stdout()

                script:
                ''
            }

            workflow {
                hello()
                hello.out
            }
            ''',
            'Using the `.out` property to access process/workflow outputs is not supported with static typing -- assign the output to a variable instead'
        )
    }

    def 'should not allow a void call result to be assigned to a variable' () {
        expect:
        check(
            '''\
            def foo() {
                println 'hello'
            }

            workflow {
                foo()
                def x = foo()
            }
            ''',
            'Cannot assign the result of a call that does not return a value'
        )

        and:
        check(
            '''\
            nextflow.enable.types = true

            process foo {
                input:
                message: String

                script:
                ''
            }

            workflow {
                def result = foo('hello')
            }
            ''',
            'Cannot assign the result of a call that does not return a value'
        )

        and:
        check(
            '''\
            nextflow.enable.types = true

            workflow foo {
                take:
                message: String

                main:
                println message
            }

            workflow {
                def result = foo( 'hello' )
            }
            ''',
            'Cannot assign the result of a call that does not return a value'
        )
    }

    def 'should resolve record type' () {
        when:
        def exp = parseExpression(
            '''\
            record(id: '1', count: 42)
            '''
        )
        def type = getType(exp)
        then:
        Types.getName(type) == 'Record {\n    id: String\n    count: Integer\n}'
        type.getField('id') != null
        type.getField('count') != null

        when:
        exp = parseExpression(
            '''\
            record(id: '1', count: 42).id
            '''
        )
        type = getType(exp)
        then:
        Types.getName(type) == 'String'
    }

    def 'should report error for record type constructor' () {
        expect:
        check(
            '''\
            record Sample {
                id: String
                reads: List<Path>
            }

            workflow {
                new Sample()
            }
            ''',
            'Record type Sample cannot be used as a constructor -- use `record()` instead'
        )
    }

    def 'should resolve record sum' () {
        when:
        def exp = parseExpression(
            '''\
            record(id: '1') + record(count: 42)
            '''
        )
        def type = getType(exp)
        then:
        Types.getName(type) == 'Record {\n    id: String\n    count: Integer\n}'
        type.getField('id') != null
        type.getField('count') != null
    }

    def 'should resolve tuple type' () {
        when:
        def exp = parseExpression(
            '''\
            tuple('hello', 42, true)
            '''
        )
        def type = getType(exp)
        then:
        Types.getName(type) == 'Tuple<String, Integer, Boolean>'

        when:
        exp = parseExpression(
            '''\
            tuple('hello', 42, true)[1]
            '''
        )
        type = getType(exp)
        then:
        Types.getName(type) == 'Integer'
    }

    def 'should resolve tuple assignment' () {
        when:
        def exp = parseExpression(
            '''\
            (x, y, z) = tuple('hello', 42, true)
            tuple(x, y, z)
            '''
        )
        def type = getType(exp)
        then:
        Types.getName(type) == 'Tuple<String, Integer, Boolean>'
    }

    def 'should allow tuple destructuring in closure' () {
        expect:
        check(
            '''\
            channel.value( tuple(1, 2, 3) ).map { x, y, z -> y }
            ''',
            null
        )
    }

    @Unroll
    def 'should check a spread-dot expression' () {
        expect:
        check(SOURCE, ERROR)

        where:
        SOURCE                          | ERROR
        "'*.txt'*.dirName"              | "Spread-dot is only supported for Iterable types"
        "files('*.txt')*.dirName"       | "Unrecognized property `dirName` for element type Path"
        "files('*.txt')*.name"          | null
        "'*.txt'*.getDirName()"         | "Spread-dot is only supported for Iterable types"
        "files('*.txt')*.getDirName()"  | "Unrecognized method `getDirName` for element type Path"
        "files('*.txt')*.toUriString()" | null
    }

    def 'should resolve a `combine` operation' () {
        when:
        def exp = parseExpression(
            '''\
            left  = channel.of( 42 )
            right = channel.of( 'hello' )
            left.combine(right)
            '''
        )
        def type = getType(exp)
        then:
        Types.getName(type) == 'Channel<Tuple<Integer, String>>'

        when:
        exp = parseExpression(
            '''\
            left  = channel.of( tuple(1, 2) )
            right = channel.of( tuple('red', 'blue') )
            left.combine(right)
            '''
        )
        type = getType(exp)
        then:
        Types.getName(type) == 'Channel<Tuple<Integer, Integer, String, String>>'

        when:
        exp = parseExpression(
            '''\
            samples = channel.of( record(id: 1, fastq: file('1.fq')) )
            index = channel.value( file('index.fa') )
            samples.combine( single_end: true, index: index )
            '''
        )
        type = getType(exp)
        then:
        Types.getName(type) == 'Channel<Record {\n    id: Integer\n    fastq: Path\n    single_end: Boolean\n    index: Path\n}>'
    }

    def 'should resolve a `combine` operation with a raw tuple' () {
        when:
        def exp = parseExpression(
            '''\
            def pair: Tuple = tuple(1, 'a')
            left  = channel.of( 42 )
            right = channel.of( pair )
            left.combine(right)
            '''
        )
        then:
        Types.getName(getType(exp)) == 'Channel<Tuple>'
    }

    def 'should resolve a `combine` operation on a dataflow value' () {
        when:
        def exp = parseExpression(
            '''\
            left  = channel.value( 42 )
            right = channel.value( 'hello' )
            left.combine(right)
            '''
        )
        def type = getType(exp)
        then:
        Types.getName(type) == 'Value<Tuple<Integer, String>>'

        when:
        exp = parseExpression(
            '''\
            sample = channel.value( record(id: 1, fastq: file('1.fq')) )
            index = channel.value( file('index.fa') )
            sample.combine( single_end: true, index: index )
            '''
        )
        type = getType(exp)
        then:
        Types.getName(type) == 'Value<Record {\n    id: Integer\n    fastq: Path\n    single_end: Boolean\n    index: Path\n}>'
    }

    def 'should report an error for a `combine` named argument that is a channel' () {
        when:
        def errors = getErrors(
            '''\
            samples = channel.of( record(id: 1, fastq: file('1.fq')) )
            index = channel.of( file('index.fa') )
            samples.combine( index: index )
            '''
        )
        then:
        errors.size() == 1
        errors[0].getOriginalMessage() == 'Named argument of a `combine` operation cannot be a channel -- it must be a dataflow value or plain value'
    }

    def 'should resolve a `groupBy` operation' () {
        when:
        def exp = parseExpression(
            '''\
            left = channel.of( tuple(42, 'hello'), tuple(42, 'goodbye') )
            left.groupBy()
            '''
        )
        def type = getType(exp)
        then:
        Types.getName(type) == 'Channel<Tuple<Integer, Bag<String>>>'

        when:
        exp = parseExpression(
            '''\
            left = channel.of( tuple(42, 2, 'hello'), tuple(42, 2, 'goodbye') )
            left.groupBy()
            '''
        )
        type = getType(exp)
        then:
        Types.getName(type) == 'Channel<Tuple<Integer, Bag<String>>>'
    }

    def 'should resolve a `join` operation' () {
        when:
        def exp = parseExpression(
            '''\
            left  = channel.of( record(id: 42, name: 'hello') )
            right = channel.of( record(id: 42, alive: true) )
            left.join(right, by: 'id')
            '''
        )
        def type = getType(exp)
        then:
        Types.getName(type) == 'Channel<Record {\n    id: Integer\n    name: String\n    alive: Boolean\n}>'
    }

    def 'should report an error for a `join` field that is not in both records' () {
        when:
        def errors = getErrors(
            '''\
            left  = channel.of( record(id: 42, name: 'hello') )
            right = channel.of( record(id: 42, alive: true) )
            left.join(right, by: 'sample_id')
            '''
        )
        then:
        errors.size() == 2
        errors[0].getOriginalMessage() == 'Join field `sample_id` is not present in left-hand side'
        errors[1].getOriginalMessage() == 'Join field `sample_id` is not present in right-hand side'
    }

}
