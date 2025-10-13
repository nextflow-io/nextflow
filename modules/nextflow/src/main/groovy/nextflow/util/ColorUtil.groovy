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

package nextflow.util

import groovy.transform.CompileStatic
import nextflow.SysEnv
import org.fusesource.jansi.Ansi
import static org.fusesource.jansi.Ansi.*

/**
 * Utility class for ANSI color formatting
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@CompileStatic
class ColorUtil {

    /**
     * Check if ANSI colors should be enabled based on Nextflow conventions
     */
    static boolean isAnsiEnabled() {
        // Check NO_COLOR environment variable (https://no-color.org/)
        if (SysEnv.get('NO_COLOR')) {
            return false
        }

        // Check NXF_ANSI_LOG environment variable
        final env = SysEnv.get('NXF_ANSI_LOG')
        if (env) {
            try {
                return Boolean.parseBoolean(env)
            } catch (Exception e) {
                // Invalid value, fall through to default
            }
        }

        // Default to enabled if terminal supports it
        return Ansi.isEnabled()
    }

    /**
     * Print colored text using format keywords and reset all styles when done
     * Format: colorize(text, format) - e.g. colorize('Hello', 'red bold'), colorize('World', 'cyan dim')
     */
    static void printColored(String text, String format = '') {
        println colorize(text, format, true) // Always do full reset for printColored
    }

    /**
     * Format text with color using format keywords
     * Format: colorize(text, format, fullReset) - e.g. colorize('Hello', 'red bold'), colorize('World', 'cyan on yellow')
     * Supported colors: black, red, green, yellow, blue, magenta, cyan, white, default
     * Supported background: on black, on red, on green, on yellow, on blue, on magenta, on cyan, on white, on default
     * Supported formatting: bold, dim
     * @param fullReset if true, resets all styles; if false, only resets styles that were applied
     */
    static String colorize(String text, String format = '', boolean fullReset = false) {
        if (!isAnsiEnabled() || !text) {
            return text ?: ''
        }

        if (!format) {
            return text
        }

        def parts = format.split(' ')
        def foregroundColor = null
        def backgroundColor = null
        def bold = false
        def dim = false
        def onNext = false

        // Parse format keywords
        parts.each { part ->
            def partLower = part.toLowerCase()

            if (onNext) {
                // Previous word was 'on', so this is a background color
                backgroundColor = parseColor(partLower)
                onNext = false
            } else if (partLower == 'on') {
                onNext = true
            } else {
                switch (partLower) {
                    case 'black':
                    case 'red':
                    case 'green':
                    case 'yellow':
                    case 'blue':
                    case 'magenta':
                    case 'cyan':
                    case 'white':
                    case 'default':
                        foregroundColor = parseColor(partLower)
                        break
                    case 'bold':
                        bold = true
                        break
                    case 'dim':
                        dim = true
                        break
                }
            }
        }

        def fmt = ansi()
        if (foregroundColor) fmt = fmt.fg(foregroundColor)
        if (backgroundColor) fmt = fmt.bg(backgroundColor)
        if (bold) fmt = fmt.bold()
        if (dim) fmt = fmt.a(Attribute.INTENSITY_FAINT)
        fmt = fmt.a(text)

        if (fullReset) {
            // Reset all styles
            fmt = fmt.a(Attribute.RESET)
        } else {
            // Reset only what was set
            if (bold) fmt = fmt.boldOff()
            if (dim) fmt = fmt.a(Attribute.INTENSITY_BOLD) // Reset dim
            if (foregroundColor) fmt = fmt.fg(Color.DEFAULT)
            if (backgroundColor) fmt = fmt.bg(Color.DEFAULT)
        }

        return fmt.toString()
    }

    /**
     * Parse color name to Jansi Color enum
     */
    private static Color parseColor(String colorName) {
        switch (colorName) {
            case 'black': return Color.BLACK
            case 'red': return Color.RED
            case 'green': return Color.GREEN
            case 'yellow': return Color.YELLOW
            case 'blue': return Color.BLUE
            case 'magenta': return Color.MAGENTA
            case 'cyan': return Color.CYAN
            case 'white': return Color.WHITE
            case 'default': return Color.DEFAULT
            default: return null
        }
    }

}
