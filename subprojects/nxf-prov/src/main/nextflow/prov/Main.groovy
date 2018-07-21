package nextflow.prov

/**
 * RO-Prov command entry point
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class Main {

    static void main(String... args) {
        CmdProv cmdProv = new CmdProv()
        cmdProv.run()

        println "Hello prov!"
    }
}
