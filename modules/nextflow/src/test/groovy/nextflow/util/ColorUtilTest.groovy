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

package nextflow.util

import nextflow.SysEnv
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture

/**
 * Test ColorUtil functionality
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class ColorUtilTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    def setup() {
        // No special setup needed
    }

    def 'should disable ANSI when NO_COLOR is set'() {
        given:
        SysEnv.push([NO_COLOR: '1'])

        when:
        def result = ColorUtil.isAnsiEnabled()

        then:
        !result

        cleanup:
        SysEnv.pop()
    }

    def 'should disable ANSI when NXF_ANSI_LOG is false'() {
        given:
        SysEnv.push([NXF_ANSI_LOG: 'false'])

        when:
        def result = ColorUtil.isAnsiEnabled()

        then:
        !result

        cleanup:
        SysEnv.pop()
    }

    def 'should enable ANSI when NXF_ANSI_LOG is true'() {
        given:
        SysEnv.push([NXF_ANSI_LOG: 'true'])

        when:
        def result = ColorUtil.isAnsiEnabled()

        then:
        result

        cleanup:
        SysEnv.pop()
    }

    def 'should handle invalid NXF_ANSI_LOG value'() {
        given:
        SysEnv.push([NXF_ANSI_LOG: 'invalid'])

        when:
        def result = ColorUtil.isAnsiEnabled()

        then:
        // Should not crash and return a boolean
        result instanceof Boolean

        cleanup:
        SysEnv.pop()
    }

    def 'should colorize text with basic colors'() {
        when:
        def result = ColorUtil.colorize('test', 'red')

        then:
        result.contains('test')
        if (ColorUtil.isAnsiEnabled()) {
            assert result.contains('\u001B[31m')  // Red color code
        }
    }

    def 'should colorize text with background colors'() {
        when:
        def result = ColorUtil.colorize('test', 'red on yellow')

        then:
        result.contains('test')
        if (ColorUtil.isAnsiEnabled()) {
            assert result.contains('\u001B[31;43m')  // Red foreground, yellow background
        }
    }

    def 'should colorize text with formatting'() {
        when:
        def result = ColorUtil.colorize('test', 'bold dim cyan')

        then:
        result.contains('test')
        if (ColorUtil.isAnsiEnabled()) {
            assert result.contains('\u001B[36;1;2m')  // Cyan color, bold, dim
        }
    }

    def 'should handle empty or null text'() {
        expect:
        ColorUtil.colorize(null, 'red') == ''
        ColorUtil.colorize('', 'red') == ''
    }

    def 'should handle empty format'() {
        when:
        def result = ColorUtil.colorize('test', '')

        then:
        result == 'test'
    }

    def 'should handle null format'() {
        when:
        def result = ColorUtil.colorize('test', null)

        then:
        result == 'test'
    }

    def 'should reset styles with fullReset=true'() {
        when:
        def result = ColorUtil.colorize('test', 'red bold', true)

        then:
        if (ColorUtil.isAnsiEnabled()) {
            assert result.endsWith('\u001B[m')  // Full reset code
        }
    }

    def 'should reset only applied styles with fullReset=false'() {
        when:
        def result = ColorUtil.colorize('test', 'red bold', false)

        then:
        if (ColorUtil.isAnsiEnabled()) {
            assert !result.endsWith('\u001B[m')  // Should not end with full reset
            assert result.endsWith('\u001B[22;39m')  // Bold off, default foreground
        }
    }

    def 'should print colored text'() {
        when:
        ColorUtil.printColored('test message', 'green')

        then:
        def output = capture.toString()
        output.contains('test message')
        if (ColorUtil.isAnsiEnabled()) {
            assert output.contains('\u001B[32m')  // Green color
            assert output.contains('\u001B[m')   // Reset at end
        }
    }

    def 'should handle all supported colors'() {
        given:
        def colors = ['black', 'red', 'green', 'yellow', 'blue', 'magenta', 'cyan', 'white', 'default']

        expect:
        colors.each { color ->
            def result = ColorUtil.colorize('test', color)
            assert result.contains('test')
        }
    }

    def 'should handle complex format combinations'() {
        when:
        def result = ColorUtil.colorize('complex', 'bold red on blue dim')

        then:
        result.contains('complex')
        if (ColorUtil.isAnsiEnabled()) {
            assert result.contains('\u001B[31;44;1;2m')  // Red foreground, blue background, bold, dim
        }
    }

    def 'should be case insensitive'() {
        when:
        def result1 = ColorUtil.colorize('test', 'RED BOLD')
        def result2 = ColorUtil.colorize('test', 'red bold')

        then:
        if (ColorUtil.isAnsiEnabled()) {
            assert result1 == result2
        } else {
            assert result1 == 'test'
            assert result2 == 'test'
        }
    }

    def 'should return plain text when ANSI disabled'() {
        given:
        SysEnv.push([NO_COLOR: '1'])

        when:
        def result = ColorUtil.colorize('plain text', 'red bold on yellow')

        then:
        result == 'plain text'

        cleanup:
        SysEnv.pop()
    }
}
