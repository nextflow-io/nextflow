package nextflow.extension

import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Session

interface LinChannelEx {
    void viewLineage(Session session, DataflowWriteChannel channel, URI uri)

    void queryLineage(Session session, DataflowWriteChannel channel, String query)
}