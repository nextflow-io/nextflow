/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.ui.console

import javax.swing.*
import javax.swing.filechooser.FileFilter
import java.nio.file.Path
import java.nio.file.Paths

import groovy.console.ui.Console
import groovy.console.ui.OutputTransforms
import groovy.transform.ThreadInterrupt
import groovy.util.logging.Slf4j
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.cli.CliOptions
import nextflow.cli.CmdInfo
import nextflow.cli.CmdRun
import nextflow.config.ConfigBuilder
import nextflow.script.ScriptBinding
import nextflow.script.ScriptFile
import nextflow.script.ScriptParser
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
/**
 * Implement a REPL console for Nextflow DSL based on Groovy console
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Nextflow extends Console {

    static public final TITLE = 'Nextflow REPL console'

    static {

        Console.groovyFileFilter = new NextflowFileFilter()
        Console.frameConsoleDelegates = [
                rootContainerDelegate:{
                    frame(
                            title: TITLE,
                            //location: [100,100], // in groovy 2.0 use platform default location
                            iconImage: imageIcon('/nextflow-icon.png').image,
                            defaultCloseOperation: JFrame.DISPOSE_ON_CLOSE,
                    ) {
                        try {
                            current.locationByPlatform = true
                        } catch (Exception e) {
                            current.location = [100, 100] // for 1.4 compatibility
                        }
                        containingWindows += current
                    }
                },
                menuBarDelegate: {arg-> current.JMenuBar = build(arg)}
        ];
    }

    private Map scriptConfig

    Nextflow(ClassLoader loader) {
        super(loader, new ScriptBinding())
        this.scriptConfig = createScriptConfig()
        NextflowMeta.instance.enableDsl2()
    }

    protected Map createScriptConfig() {
        final script = scriptFile as Path
        final base = script ? script.parent : Paths.get('.')

        // create the config object
        return new ConfigBuilder()
                    .setOptions( new CliOptions() )
                    .setBaseDir(base)
                    .setCmdRun( new CmdRun() )
                    .build()
    }

    /**
     * Create a new Groovy Shell to interpret Nextflow scripts
     *
     * @param parent
     * @param binding
     */
    @Override
    void newScript(ClassLoader parent, Binding binding) {
        assert parent

        def parser = new ScriptParser(parent)
        config = parser.getConfig()

        if (threadInterrupt)
            config.addCompilationCustomizers(new ASTTransformationCustomizer(ThreadInterrupt))

        parser.setBinding((ScriptBinding)binding)
        shell = parser.getInterpreter()
    }

    @Override
    void clearContext(EventObject evt = null) {
        final binding = new ScriptBinding()
        newScript(null, binding)
        // reload output transforms
        binding.variables._outputTransforms = OutputTransforms.loadOutputTransforms()
    }

    @Override
    void runScript(EventObject event = null) {
        runWith { super.runScript(event) }
    }

    @Override
    void runSelectedScript(EventObject event = null) {
        runWith { super.runSelectedScript(event) }
    }

    private runWith( Closure launcher ) {

        def script = scriptFile ? new ScriptFile((File)scriptFile) : null
        def path = scriptFile as Path
        def session = new Session(scriptConfig).init(script)
        def binding = (ScriptBinding)shell.getContext()

        binding.setSession(session)
        binding.setScriptPath(path)

        beforeExecution = {
            session.start()
        }
        afterExecution = {
            session.fireDataflowNetwork()
            session.await()
            session.destroy()
        }

        launcher.call()
    }


    /**
     * Update the window title
     */
    @Override
    void updateTitle() {
        if (frame.properties.containsKey('title')) {
            if (scriptFile != null) {
                frame.title = scriptFile.name + (dirty?' * ':'') + ' - ' + TITLE
            } else {
                frame.title = TITLE
            }
        }
    }

    /**
     * Show a customized about dialog box
     *
     * @param evt
     */
    void showAbout(EventObject evt = null) {
        def pane = swing.optionPane()
        pane.setMessage('REPL Console for evaluating Nextflow scripts\n\n' + CmdInfo.getInfo(0))
        def dialog = pane.createDialog(frame, 'About ' + TITLE)
        dialog.show()
    }


    @Override
    def finishNormal(Object result) {
        // Take down the wait/cancel dialog
        history[-1].result = result
        statusLabel.text = 'Execution complete.'

        if( !visualizeScriptResults )
            return

        if (result != null) {
            appendOutputNl('Result: ', promptStyle)
            def obj = OutputTransforms.transformResult(result, shell.context._outputTransforms)

            // multi-methods are magical!
            appendOutput(obj, resultStyle)
        } else {
            statusLabel.text = 'Execution complete. Result was null.'
        }
        bindResults()
        if (detachedOutput) {
            prepareOutputWindow()
            showOutputWindow()
        }
    }


    /**
     * Filter supported script files in open dialog
     */
    private static class NextflowFileFilter extends FileFilter {
        private static final SOURCE_EXTENSIONS = ['*.nf', '*.groovy']
        private static final SOURCE_EXT_DESC = SOURCE_EXTENSIONS.join(',')

        boolean accept(File f) {
            if (f.isDirectory()) {
                return true
            }
            SOURCE_EXTENSIONS.find {it == getExtension(f)} ? true : false
        }

        String getDescription() {
            "Nextflow scripts ($SOURCE_EXT_DESC)"
        }

        static String getExtension(f) {
            def ext = null;
            def s = f.getName()
            def i = s.lastIndexOf('.')
            if (i > 0 &&  i < s.length() - 1) {
                ext = s.substring(i).toLowerCase()
            }
            "*$ext"
        }
    }

}
