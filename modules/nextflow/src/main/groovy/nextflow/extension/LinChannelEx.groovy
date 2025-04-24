package nextflow.extension

import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Session

interface LinChannelEx {
    void queryLineage(Session session, DataflowWriteChannel channel, URI uri)

}