package nextflow.sql


import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.NF
import nextflow.Session
import nextflow.extension.CH
import nextflow.extension.ChannelExtensionPoint
/**
 * Provide a channel factory extension that allows the execution of Sql queries
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ChannelSqlExtension implements ChannelExtensionPoint {

    private Session session
    private Map config

    void init(Session session) {
        this.session = session
        this.config = (Map) session.config.dataSources
    }

    @Override
    String getScope() {
        return 'sql'
    }

    DataflowWriteChannel fromQuery(String query) {
        fromQuery([:], query)
    }

    DataflowWriteChannel fromQuery(Map opts, String query) {
        return queryToChannel(query, opts)
    }

    protected DataflowWriteChannel queryToChannel(String query, Map opts) {
        final channel = CH.create()
        if(NF.dsl2) {
            session.addIgniter {-> new QueryOpImpl().withOpts(opts).withStatement(query).withTarget(channel).perform(true) }
        }
        else {
            new QueryOpImpl()
                    .withOpts(opts)
                    .withStatement(query)
                    .withTarget(channel)
                    .perform(true)
        }
        return channel
    }

}
