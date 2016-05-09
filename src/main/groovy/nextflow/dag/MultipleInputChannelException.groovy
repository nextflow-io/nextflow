package nextflow.dag
import groovy.transform.ToString
import nextflow.dag.DAG.ChannelHandler
import nextflow.dag.DAG.Vertex
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true)
class MultipleInputChannelException extends Exception {

    String name

    ChannelHandler channel

    Vertex duplicate

    Vertex existing

    MultipleInputChannelException( String name, ChannelHandler channel, Vertex duplicate, Vertex existing) {
        super()
        this.name = name
        this.channel = channel
        this.duplicate = duplicate
        this.existing = existing
    }

    String getMessage() {
        if( !name ) {
            return 'Channels cannot be used as input in more than one process or operator'
        }

        if( duplicate.type != DAG.Type.PROCESS ) {
            return "Channel `$name` has been used as an input by more than a process or an operator"
        }

        String message = "Channel `$name` has been used twice as an input by process `${duplicate.label}`"
        if( existing.type == DAG.Type.PROCESS )
            message += " and process `${existing.label}`"
        else
            message += " and another operator"
        return message
    }

}
