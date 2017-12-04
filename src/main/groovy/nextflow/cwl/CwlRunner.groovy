package nextflow.cwl
import java.nio.file.Path

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.processor.ProcessConfig
import nextflow.script.BaseScript
import nextflow.script.TaskBody
import nextflow.script.TokenVar
import org.yaml.snakeyaml.Yaml
/**
 * Model a CWL step. This class fetch a CWL CommandLineTool descriptor,
 * parse it, and creates the corresponding {@link nextflow.processor.ProcessConfig} object
 * representing a Nextflow process configuration (directives, inputs and outputs)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class CwlRunner {

    private BaseScript owner

    String command

    ProcessConfig config

    Map<String,Object> inputs

    Binding outputs

    CwlRunner( BaseScript owner, Map inputs = [:] ) {
        this.owner = owner
        this.inputs = inputs
        this.outputs = new Binding()
        this.config = new ProcessConfig(owner)
        this.config.setBinding(outputs)
        this.config.placeholder = '$'
    }


    def parse( Map yaml ) {
        if( !yaml.cwlVersion )
            throw new IllegalArgumentException("Not a valid CWL descriptor -- Missing `cwlVersion`")
        if( yaml.class != 'CommandLineTool' )
            throw new IllegalArgumentException("Not a valid CWL descriptor -- Invalid class: ${yaml.class}")

        def cmd = new StringBuilder()
        // -- base command
        if( yaml.baseCommand )
            cmd << yaml.baseCommand
        else
            throw new IllegalArgumentException("Not a valid CWL descriptor -- Missing `baseCommand`")

        // -- command arguments
        def args = (List)yaml.arguments
        if( args )
            cmd << ' ' << args.join(' ')

        // collect inputs
        yaml.inputs.each { String name, Map value ->
            cmd << ' ' << '${' << name << '}'
            createInputParam(name, value)
        }

        // collect outputs
        yaml.outputs.each { String name, Map value ->
            createOutputParam(name, value)
        }

        this.command = cmd.toString()
    }

    def parse( String yaml ) {
        parse((Map)new Yaml().load(yaml))
    }

    protected void createInputParam(String name, Map value) {
        if( !inputs.containsKey(name) )
            throw new IllegalArgumentException("CWL missing input parameter: $name")

        final ch = inputs.get(name)
        final type = value.type?.toString()?.toLowerCase()
        if( type.startsWith('file') ) {
            config._in_file(new TokenVar(name)).from(ch)
        }
        else if ( type.startsWith('string') ) {
            config._in_val(new TokenVar(name)).from(ch)
        }
        else {
            throw new IllegalArgumentException("CWL unknown input type: $type")
        }
    }

    protected void createOutputParam(String name, Map value) {
        final type = value.type?.toString()?.toLowerCase()

        if( type.startsWith('file') ) {
            def bind = value.outputBinding
            def glob = bind instanceof Map ? bind.glob : null
            if( !glob ) throw new IllegalArgumentException("CWL missing output file glob")
            config._out_file(glob).into(new TokenVar(name))
        }
//        else if( type.startsWith('string') ) {
//
//        }
        else {
            throw new IllegalArgumentException("Unknown CWL output type: $type")
        }
    }


    Map<String,DataflowWriteChannel> run(Path path) {
        // parse the CWL command
        def source = path.text
        parse(source)

        def body = new TaskBody(new CwlCommand(command), source, 'shell')
        def proc = owner.getProcessFactory().newTaskProcessor('name', config, body)
        proc.run()
        def outChannels = outputs.getVariables()
        log.debug "CWL task output channels: $outChannels"
        return outChannels
    }

}