/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import groovy.transform.Field
import groovy.transform.ThreadInterrupt
import groovy.ui.Console
import groovy.ui.OutputTransforms
import groovy.util.logging.Slf4j
import nextflow.Channel
import nextflow.Session
import nextflow.ast.NextflowDSL
import nextflow.ast.NextflowXform
import nextflow.cli.CliOptions
import nextflow.cli.CmdInfo
import nextflow.cli.CmdRun
import nextflow.script.BaseScript
import nextflow.config.ConfigBuilder
import nextflow.script.ScriptBinding
import nextflow.util.Duration
import nextflow.util.LoggerHelper
import nextflow.util.MemoryUnit
import org.apache.commons.lang.StringUtils
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.control.customizers.ImportCustomizer
import org.codehaus.groovy.runtime.StackTraceUtils
/**
 * Implement a REPL console for Nextflow DSL based on Groovy console
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Nextflow extends Console {

    static final TITLE = 'Nextflow REPL console'

    static {

        Console.groovyFileFilter = new NextflowFileFilter()
        Console.frameConsoleDelegates = [
                rootContainerDelegate:{
                    frame(
                            title: TITLE,
                            //location: [100,100], // in groovy 2.0 use platform default location
                            iconImage: imageIcon('/nextflow_icon_48x48.png').image,
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

    private ClassLoader loader

    Nextflow(ClassLoader loader) {
        super(loader, createBinding())
        this.loader = loader
    }

    static protected Binding createBinding() {

        // create the config object
        def config = new ConfigBuilder()
                        .setOptions( new CliOptions() )
                        .setBaseDir(Paths.get('.'))
                        .setCmdRun( new CmdRun() )
                        .build()

        return new ScriptBinding(config)
    }


    protected CompilerConfiguration createCompilerConfig() {

        // define the imports
        def importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( StringUtils.name, Field.name )
        importCustomizer.addImports( Path.name )
        importCustomizer.addImports( Channel.name )
        importCustomizer.addImports( Duration.name )
        importCustomizer.addImports( MemoryUnit.name )
        importCustomizer.addStaticStars( nextflow.Nextflow.name )

        def config = new CompilerConfiguration()
        config.addCompilationCustomizers( importCustomizer )
        config.scriptBaseClass = BaseScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowXform))

        return config
    }

    /**
     * Create a new Groovy Shell to interpret Nextflow scripts
     *
     * @param parent
     * @param binding
     */
    @Override
    void newScript(ClassLoader parent, Binding binding) {

        config = createCompilerConfig()

        // from super 'newScript' implementation
        if (threadInterrupt) config.addCompilationCustomizers(new ASTTransformationCustomizer(ThreadInterrupt))

        // run and wait for termination
        shell = new GroovyShell(loader, binding, config)
    }

    @Override
    void clearContext(EventObject evt = null) {
        def binding = createBinding()
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

        def binding = shell.context as ScriptBinding
        def session = new Session(binding)
        session.init(scriptFile as Path)

        binding.setVariable( 'baseDir', session.baseDir )
        binding.setVariable( 'workDir', session.workDir )

        beforeExecution = { session.start() }
        afterExecution = { session.await(); session.destroy() }

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
     * Nextflow REPL entry point
     *
     * @param args
     */
    static void main(String... args) {
        CliOptions opts = new CliOptions()
        opts.logFile = '.nextflow-console.log'
        new LoggerHelper(opts).setup()

        if (args.length == 2 && args[1] == '--help') {
            println 'usage: nextflow console [filename]'
            return
        }

        // full stack trace should not be logged to the output window - GROOVY-4663
        java.util.logging.Logger.getLogger(StackTraceUtils.STACK_LOG_NAME).useParentHandlers = false

        //when starting via main set the look and feel to system
        UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName())

        def console = new Nextflow(Nextflow.class.classLoader?.getRootLoader())
        console.useScriptClassLoaderForScriptExecution = true
        console.run()
        if (args.length == 2)
            try {
                console.loadScriptFile(args[1] as File)
            }
            catch( IOException e ) {
                log.warn("Can't open script file: ${args[1]}" )
            }

    }

    /**
     * Filter supported script files in open dialog
     */
    private static class NextflowFileFilter extends FileFilter {
        private static final SOURCE_EXTENSIONS = ['*.nf', '*.groovy']
        private static final SOURCE_EXT_DESC = SOURCE_EXTENSIONS.join(',')

        public boolean accept(File f) {
            if (f.isDirectory()) {
                return true
            }
            SOURCE_EXTENSIONS.find {it == getExtension(f)} ? true : false
        }

        public String getDescription() {
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
