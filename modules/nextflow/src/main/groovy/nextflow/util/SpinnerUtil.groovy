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
import org.fusesource.jansi.Ansi
import org.fusesource.jansi.AnsiConsole

import java.util.concurrent.atomic.AtomicBoolean

import static org.fusesource.jansi.Ansi.Erase

/**
 * Spinner utility for showing animated progress indicators
 *
 * Usage:
 * <pre>
 * def spinner = new SpinnerUtil("Loading...")
 * spinner.start()
 * try {
 *     // Do work
 *     spinner.updateMessage("Still loading...")
 * } finally {
 *     spinner.stop()
 * }
 * </pre>
 *
 * Or use the withSpinner method for automatic cleanup:
 * <pre>
 * SpinnerUtil.withSpinner("Loading...") { spinner ->
 *     // Do work
 *     spinner.updateMessage("Still loading...")
 * }
 * </pre>
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@CompileStatic
class SpinnerUtil {

    private static final String[] SPINNER_CHARS_WAITING = ['⠋', '⠉', '⠙', '⠘', '⠐', '⠘', '⠙', '⠉', '⠋', '⠃', '⠂', '⠃'] as String[]
    private static final String[] SPINNER_CHARS_RUNNING = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏'] as String[]
    private static final String[] SPINNER_CHARS_SUCCEEDED = ['⠀','⠄','⠤','⠦','⠶','⠷','⠿','⠻','⠛','⠙','⠉','⠈'] as String[]
    private static final String[] SPINNER_CHARS_FAILED = ['⠈', '⠉', '⠙', '⠛', '⠻', '⠿', '⠷', '⠶', '⠦', '⠤', '⠄', '⠀'] as String[]
    private static final int SPINNER_UPDATE_MS = 100 // Update spinner every 100ms for smooth animation

    private final AtomicBoolean shouldStop = new AtomicBoolean(false)
    private String message
    private String colorName
    private String[] spinnerChars
    private Thread spinnerThread
    private final boolean ansiEnabled
    private volatile int spinnerIndex = 0

    /**
     * Create a new spinner with the given message
     *
     * @param message The message to display next to the spinner
     * @param waiting If true, use waiting animation; if false, use running animation (default)
     */
    SpinnerUtil(String message, boolean waiting = false) {
        this.message = message
        this.colorName = 'cyan'
        this.spinnerChars = waiting ? SPINNER_CHARS_WAITING : SPINNER_CHARS_RUNNING
        this.ansiEnabled = ColorUtil.isAnsiEnabled()
    }

    /**
     * Start the spinner animation in a background thread
     */
    void start() {
        if (!ansiEnabled || spinnerThread != null) {
            return
        }

        shouldStop.set(false)
        spinnerThread = new Thread({
            try {
                while (!shouldStop.get()) {
                    def fmt = Ansi.ansi()
                    fmt.a("\r").eraseLine(Erase.ALL)
                    fmt.fg(parseColor(colorName)).a(spinnerChars[spinnerIndex]).reset()
                    fmt.a(" ").a(message)
                    AnsiConsole.out.print(fmt.toString())
                    AnsiConsole.out.flush()
                    spinnerIndex = (spinnerIndex + 1) % spinnerChars.length
                    Thread.sleep(SPINNER_UPDATE_MS)
                }
            } catch (InterruptedException e) {
                // Thread interrupted, exit gracefully
            } catch (Exception e) {
                // Ignore other errors to prevent spinner from breaking the main application
            }
        })
        spinnerThread.daemon = true
        spinnerThread.start()
    }

    /**
     * Update the message displayed next to the spinner
     *
     * @param newMessage The new message to display
     */
    void updateMessage(String newMessage) {
        this.message = newMessage
    }

    /**
     * Update the message and color of the spinner
     *
     * @param newMessage The new message to display
     * @param newColorName The new color name for the spinner character (e.g., 'yellow', 'blue', 'red')
     */
    void updateMessage(String newMessage, String newColorName) {
        this.message = newMessage
        this.colorName = newColorName
    }

    /**
     * Update the message, color, and mode of the spinner
     *
     * @param newMessage The new message to display
     * @param newColorName The new color name for the spinner character (e.g., 'yellow', 'blue', 'red')
     * @param waiting If true, use waiting animation; if false, use running animation
     */
    void updateMessage(String newMessage, String newColorName, boolean waiting) {
        this.message = newMessage
        this.colorName = newColorName
        this.spinnerChars = waiting ? SPINNER_CHARS_WAITING : SPINNER_CHARS_RUNNING
        this.spinnerIndex = 0 // Reset animation when changing mode
    }

    /**
     * Update the message, color, and mode of the spinner with string mode
     *
     * @param newMessage The new message to display
     * @param newColorName The new color name for the spinner character (e.g., 'yellow', 'blue', 'red')
     * @param mode The spinner mode: 'waiting', 'running', 'succeeded', or 'failed'
     */
    void updateMessage(String newMessage, String newColorName, String mode) {
        this.message = newMessage
        this.colorName = newColorName

        switch (mode) {
            case 'waiting':
                this.spinnerChars = SPINNER_CHARS_WAITING
                break
            case 'succeeded':
                this.spinnerChars = SPINNER_CHARS_SUCCEEDED
                break
            case 'failed':
                this.spinnerChars = SPINNER_CHARS_FAILED
                break
            case 'running':
            default:
                this.spinnerChars = SPINNER_CHARS_RUNNING
                break
        }

        this.spinnerIndex = 0 // Reset animation when changing mode
    }

    /**
     * Parse color name to Jansi Color enum
     */
    private static Ansi.Color parseColor(String colorName) {
        switch (colorName) {
            case 'black': return Ansi.Color.BLACK
            case 'red': return Ansi.Color.RED
            case 'green': return Ansi.Color.GREEN
            case 'yellow': return Ansi.Color.YELLOW
            case 'blue': return Ansi.Color.BLUE
            case 'magenta': return Ansi.Color.MAGENTA
            case 'cyan': return Ansi.Color.CYAN
            case 'white': return Ansi.Color.WHITE
            case 'default': return Ansi.Color.DEFAULT
            default: return Ansi.Color.CYAN
        }
    }

    /**
     * Stop the spinner and clear the line
     *
     * @param clearLine If true, clears the spinner line completely. If false, leaves the last message visible
     */
    void stop(boolean clearLine = true) {
        shouldStop.set(true)

        if (spinnerThread != null && spinnerThread.isAlive()) {
            try {
                spinnerThread.join(200) // Wait up to 200ms for spinner to stop
            } catch (InterruptedException e) {
                // Ignore
            }
        }

        if (ansiEnabled && clearLine) {
            def fmt = Ansi.ansi()
            fmt.a("\r").eraseLine(Erase.ALL)
            AnsiConsole.out.print(fmt.toString())
            AnsiConsole.out.flush()
        }

        spinnerThread = null
    }

    /**
     * Stop the spinner and replace it with a final message
     *
     * @param finalMessage The final message to display (will be printed on a new line)
     */
    void stopWithMessage(String finalMessage) {
        shouldStop.set(true)

        if (spinnerThread != null && spinnerThread.isAlive()) {
            try {
                spinnerThread.join(200)
            } catch (InterruptedException e) {
                // Ignore
            }
        }

        if (ansiEnabled) {
            def fmt = Ansi.ansi()
            fmt.a("\r").eraseLine(Erase.ALL)
            fmt.a(finalMessage).a("\n")
            AnsiConsole.out.print(fmt.toString())
            AnsiConsole.out.flush()
        } else {
            println finalMessage
        }

        spinnerThread = null
    }

    /**
     * Check if the spinner is currently running
     */
    boolean isRunning() {
        return spinnerThread != null && spinnerThread.isAlive() && !shouldStop.get()
    }

    /**
     * Execute a closure with a spinner, automatically cleaning up when done
     *
     * @param message The message to display
     * @param closure The code to execute while the spinner is running
     */
    static void withSpinner(String message, Closure closure) {
        def spinner = new SpinnerUtil(message)
        spinner.start()
        try {
            closure.call(spinner)
        } finally {
            spinner.stop()
        }
    }
}
