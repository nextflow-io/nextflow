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

package nextflow.file

import java.nio.file.FileSystems
import java.nio.file.Paths
import java.util.regex.Pattern
import java.util.regex.PatternSyntaxException

import spock.lang.Specification
import spock.lang.Unroll

/**
 * Differential tests of {@link Globs} against the default file system glob matcher
 */
class GlobsTest extends Specification {

    static final List<String> PATHS = [
            'a', 'b', 'c', 'z', 'A', 'Z', '0', 'x', 'ab', 'ac', 'bc', 'ad', 'bd', 'abc', 'aaa', 'ax', 'xa', 'axb', 'a*b',
            '/', '/a', '/a/b', 'a/b', 'a/b/c', 'b/a', 'x/a', 'a/x', 'x/y/a', 'a/b.txt', '/x/a/b.txt',
            '.', '..', '.x', 'a.x', 'a.b', 'a.txt', 'b.txt', 'a.txt.gz', 'dir/a.txt', 'dir/sub/a.txt', '.hidden', 'dir/.hidden',
            '-', '!', '^', '&', '%', '+', ',', '\\', '`', '_', '*', '?', ' ', ':', '$', '|', '(', ')', '#',
            '[', ']', '{', '}', '[]', '[a', 'a]', ']a', '{a', 'a}', '}a', 'a{b', 'a{b,c}', 'a,b', 'b,c',
            '-a', 'a-', '+]', '-]', 'a+', '(a)', 'a|b', 'a$', '^a', '\\a', 'a\\', 'a\\b', '\\\\', '\\]',
            'é', '日本', '日', '😀', '😁', 'a😀', '￿',
            '\n', '\r', '\u0085', ' ', 'a\nb', 'x/\n', '\n/x', '\na',
    ]

    static final List<String> GLOBS = [
            // plain and wildcards
            'a', 'abc', 'a/b', '/a', '*', '**', '***', '****', '?', '??', '*a*', '*.txt', '**.txt', '**/*.txt', 'a/**', '**/a',
            'a*/b', '*/', 'dir/*', 'dir/**', '?**', '**?', '***a', 'a***', '/**', '*.*', '.*', '*\\**',
            // regex metacharacters are literal
            'a.b', 'a+', '(a)', 'a|b', 'a$', '^a', '$', '|', '#', 'a,b', '}', '}a', '{a,b}}', ']', 'a]',
            // escapes
            '\\*', '\\?', '\\[', '\\]', '\\{', '\\}', '\\,', '\\a', '\\\\', 'a\\\\b', '\\/', '\\😀', '\\é',
            // character classes
            '[a]', '[abc]', '[a-c]', '[!a]', '[!a-c]', '[^a]', '[^a-c]', '[^]', '[^^]', '[!^]', '[!^-]', '[^!]', '[^!a]', '[^!-a]', '[!^-a]',
            '[a-]', '[-a]', '[-]', '[!-]', '[!-a]', '[ab-]', '[!a-]', '[a-a]', '[!a-a]', '[+-.]', '[+--]', '[+--a]', '[.-0]', '[.-/]', '[!.-0]',
            '[a-cA-C]', '[a-c]x', '[!-a-c]', '[a-c!]', '[!!]', '[!!-a]', '[a!]', '[[]', '[a[]', '[[-a]', '[[:alpha:]]',
            '[a&&b]', '[&&]', '[a&]', '[&&a]', '[a&&&b]', '[%-&]', '[%-&&]', '[%-&&a]', '[!%-&&a]', '[&-a]', '[a-z&&[^b]]', '[A-Z&&b]',
            '[\\]', '[\\\\]', '[\\a]', '[a\\]]', '[\\]]', '[a\\-c]', '[\\-]', '[\\-a]', '[A-^]', '[A-^a]', '[A-^^]', '[A-`]', '[ -!]',
            '[*]', '[?]', '[{]', '[}]', '[,]', '[a,{]', '[a{]', '[a-c]]', '[+-]]', '[a-]]', '[-]]', '[^]]',
            '[😀]', '[!😀]', '[a-😀]', '[日本]', '[é]',
            // groups
            '{a,b}', '{a}', '{}', '{,}', '{a,}', '{,a}', '{a,b,}', '{a,b}c', '{a,b}{c,d}', '{a,*}', '{*.x,y}', '{*/a,b}', '{**,b}',
            '{a\\,b}', '{[a,b],c}', '{a,[{]}', '{a,[,]}', '{a,[}]}', '{a,b}*', '*.{txt,gz}', '**/*.{txt,x}', '{}}', '{a,?}',
    ]

    static final List<String> INVALID_GLOBS = [
            '[', '[!', '[a', '[a-', '[a\\', '[ab', '[!a', '*[', '[a-b][', '[]', '[!]', '[]]', '[!]]', '[]a]', '[!]a]', '[]-a]',
            '[z-a]', '[a-c-e]', '[a-c-]', '[a-b-]', '[a--]', '[--]', '[---]', '[!--]', '[!---]', '[--0]', '[--a]', '[^-]', '[^-a]', '[^-a-c]',
            '[A-^-]', '[a-a-]', '[!a-c-]', '[\\--a]', '[a-\\]', '[a-[]', '[😀-😂]',
            '[Z-\\]', '[Z-\\]]', '[Z-\\]a]', '[Z-\\a]', '[Z-\\\\]', '[Z-\\][a]', '[Z-[]', '[Z-[]]', '[Z-[a]]', '[Z-[][a]',
            '[/]', '[a/b]', '[!/]', '[/-0]',
            '{', '{{', '{a', '{a,b', 'a{b', '{a,b}{', '{a,{b}}', '{a\\}', '{\\}', '{a,[}', '{a,[b}]',
            '\\', 'a\\', '{a,\\',
    ]

    private static boolean jdkMatches(String glob, String path) {
        FileSystems.getDefault().getPathMatcher("glob:$glob").matches(Paths.get(path))
    }

    private static boolean globsMatches(String glob, String path) {
        Pattern.compile(Globs.toUnixRegexPattern(glob)).matcher(path).matches()
    }

    private static boolean jdkRejects(String glob) {
        try {
            FileSystems.getDefault().getPathMatcher("glob:$glob")
            return false
        }
        catch( PatternSyntaxException e ) {
            return true
        }
    }

    private static boolean globsRejects(String glob) {
        try {
            Pattern.compile(Globs.toUnixRegexPattern(glob))
            return false
        }
        catch( PatternSyntaxException e ) {
            return true
        }
    }

    def 'should use only paths which are not normalised' () {
        expect:
        PATHS.every { Paths.get(it).toString() == it }
    }

    @Unroll
    def 'should match the same paths as the JDK for glob: #glob' () {
        given:
        def mismatches = PATHS.findAll { jdkMatches(glob, it) != globsMatches(glob, it) }

        expect:
        !jdkRejects(glob)
        !globsRejects(glob)
        mismatches == []

        where:
        glob << GLOBS
    }

    @Unroll
    def 'should reject invalid glob: #glob' () {
        expect:
        jdkRejects(glob)
        globsRejects(glob)

        where:
        glob << INVALID_GLOBS
    }

    def 'should report the position of a syntax error' () {
        when:
        Globs.toUnixRegexPattern('ab[c')
        then:
        def e = thrown(PatternSyntaxException)
        e.pattern == 'ab[c'
        e.index == 2
    }

    def 'should match the JDK for random globs and paths' () {
        given:
        final random = new Random(20261028)
        final globChars = ['*', '?', '[', ']', '{', '}', ',', '!', '-', '^', '\\', '/', '&', 'a', 'b', '.', '😀'] as List<String>
        final pathChars = ['a', 'b', 'c', '.', '/', '-', '!', '^', '[', ']', '{', '}', ',', '\\', '&', '*', '😀', '\n'] as List<String>
        final paths = (1..300)
                .collect { randomString(random, pathChars, 1 + random.nextInt(6)) }
                .findAll { Paths.get(it).toString() == it }
                .unique()
        int rejected = 0
        int matched = 0
        final mismatches = []

        when:
        for( int i = 0; i < 20000; i++ ) {
            final glob = randomString(random, globChars, 1 + random.nextInt(8))
            final jdkError = jdkRejects(glob)
            if( jdkError != globsRejects(glob) ) {
                mismatches << "glob '$glob': JDK rejects=$jdkError"
                continue
            }
            if( jdkError ) {
                rejected++
                continue
            }
            final jdk = FileSystems.getDefault().getPathMatcher("glob:$glob")
            final ours = Pattern.compile(Globs.toUnixRegexPattern(glob))
            // a path resembling the glob, to make matches more likely
            final derived = glob.replace('*', 'a').replace('?', 'b').replaceAll(/[\[\]{}\\,!^]/, '')
            final candidates = (1..20).collect { paths[random.nextInt(paths.size())] }
            if( derived && Paths.get(derived).toString() == derived )
                candidates << derived
            for( String path : candidates ) {
                final expected = jdk.matches(Paths.get(path))
                if( expected )
                    matched++
                if( expected != ours.matcher(path).matches() )
                    mismatches << "glob '$glob' path '$path': JDK matches=$expected"
            }
        }

        then:
        mismatches == []
        rejected > 1000
        matched > 5000
    }

    @Unroll
    def 'should match path #PATH with glob #GLOB using the default path matcher' () {
        given:
        def matcher = FileHelper.getDefaultPathMatcher('glob:' + GLOB)

        expect:
        matcher.matches(Paths.get(PATH)) == EXPECTED

        where:
        GLOB                | PATH                  | EXPECTED
        '*.txt'             | 'foo.txt'             | true
        '*.txt'             | 'dir/foo.txt'         | false
        '**.txt'            | 'dir/foo.txt'         | true
        '**/*.txt'          | 'a/b/foo.txt'         | true
        'file_{1,2}.fq'     | 'file_2.fq'           | true
        'file_{1,2}.fq'     | 'file_3.fq'           | false
        '[!a-c]x'           | 'dx'                  | true
        '[!a-c]x'           | 'bx'                  | false
        '?x'                | '/x'                  | false
    }

    private static String randomString(Random random, List<String> alphabet, int length) {
        final result = new StringBuilder()
        for( int i = 0; i < length; i++ )
            result.append(alphabet[random.nextInt(alphabet.size())])
        return result.toString()
    }

}
