package nextflow.dag
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

class MultipleOutputChannelException extends Exception {

    String name

    DAG.ChannelHandler channel

    DAG.Vertex duplicate

    DAG.Vertex existing

    MultipleOutputChannelException( String name, DAG.ChannelHandler channel, DAG.Vertex duplicate, DAG.Vertex existing) {
        super()
        this.name = name
        this.channel = channel
        this.duplicate = duplicate
        this.existing = existing
    }

    String getMessage() {
        if( !name ) {
            return 'Channels cannot be used as output in more than one process or operator'
        }

        if( duplicate.type != DAG.Type.PROCESS ) {
            return "Channel `$name` has been used as an output by more than a process or an operator"
        }

        String message = "Channel `$name` has been used twice as an output by process `${duplicate.label}`"
        if( existing.type == DAG.Type.PROCESS )
            message += " and process `${existing.label}`"
        else
            message += " and another operator"
        return message
    }

}
