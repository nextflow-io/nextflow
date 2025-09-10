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

package nextflow.cli

import groovy.transform.CompileStatic
import nextflow.SysEnv
import org.fusesource.jansi.Ansi
import static org.fusesource.jansi.Ansi.*

/**
 * Utility class for ANSI color formatting in auth commands
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@CompileStatic
class AuthColorUtil {

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
     * Print colored text if ANSI is enabled, otherwise plain text
     */
    static void printColored(Object text, Color color = null, boolean bold = false, boolean dim = false) {
        def textStr = text?.toString() ?: ''
        if (!isAnsiEnabled() || (!color && !bold && !dim)) {
            println textStr
            return
        }
        
        def fmt = ansi()
        if (color) fmt = fmt.fg(color)
        if (bold) fmt = fmt.bold()
        if (dim) fmt = fmt.a(Attribute.INTENSITY_FAINT)
        fmt = fmt.a(textStr)
        if (dim) fmt = fmt.a(Attribute.RESET)
        if (bold) fmt = fmt.boldOff()
        if (color) fmt = fmt.fg(Color.DEFAULT)
        println fmt
    }
    
    /**
     * Format text with color inline
     */
    static String colorize(Object text, Color color = null, boolean bold = false, boolean dim = false) {
        def textStr = text?.toString() ?: ''
        if (!isAnsiEnabled() || (!color && !bold && !dim)) {
            return textStr
        }
        
        def fmt = ansi()
        if (color) fmt = fmt.fg(color)
        if (bold) fmt = fmt.bold()
        if (dim) fmt = fmt.a(Attribute.INTENSITY_FAINT)
        fmt = fmt.a(textStr)
        if (dim) fmt = fmt.a(Attribute.RESET)
        if (bold) fmt = fmt.boldOff()
        if (color) fmt = fmt.fg(Color.DEFAULT)
        return fmt.toString()
    }
}