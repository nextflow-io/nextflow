package nextflow.prov

import nextflow.cli.CliOptions
import nextflow.util.LoggerHelper

/**
 * RO-Prov command entry point
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class Main {

    static void main(String... args) {
        LoggerHelper.configureLogger(new CliOptions(logFile: '.nextflow-prov.log'))
        CmdProv cmdProv = new CmdProv()
        cmdProv.run()

    }
}
