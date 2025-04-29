package nextflow.extension

import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Session

/**
 * Interface to implement the Lineage channel factories and functions.
 * @author Jorge Ejarque <jorge.ejarque@seqera.io
 */
interface LinChannelEx {
    /**
     * Lineage metadata view.
     *
     * @param session Nextflow Session
     * @param lid Lineage Id to view
     * @return Lineage metadata object
     */
    Object viewLineage(Session session, String lid)

    /**
     * Query Lineage metadata.
     *
     * @param session Nextflow Session
     * @param channel Channel to publish the Lineage Ids matching the query params
     * @param params Query parameters
     */
    void queryLineage(Session session, DataflowWriteChannel channel, Map<String,String> params)
}