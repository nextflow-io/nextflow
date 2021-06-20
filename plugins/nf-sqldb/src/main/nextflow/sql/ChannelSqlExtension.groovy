package nextflow.sql

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.NF
import nextflow.Session
import nextflow.extension.CH
import nextflow.extension.ChannelExtensionPoint
import nextflow.sql.config.SqlConfig
/**
 * Provide a channel factory extension that allows the execution of Sql queries
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ChannelSqlExtension implements ChannelExtensionPoint {

    private Session session
    private SqlConfig config

    void init(Session session) {
        this.session = session
        this.config = new SqlConfig((Map) session.config.dataSources)
    }

    @Override
    String getScope() {
        return 'sql'
    }

    DataflowWriteChannel fromQuery(String query) {
        fromQuery(Collections.emptyMap(), query)
    }

    DataflowWriteChannel fromQuery(Map opts, String query) {
        return queryToChannel(query, opts)
    }

    protected DataflowWriteChannel queryToChannel(String query, Map opts) {
        final channel = CH.create()
        final dsName = (opts?.dataSource ?: 'default') as String
        final dataSource = config.getDatasource(dsName)
        if( dataSource==null )
            throw new IllegalArgumentException("Unknown dataSource name: $dsName")
        final handler = new QueryHandler().withDatasource(dataSource).withStatement(query).withTarget(channel)
        if(NF.dsl2) {
            session.addIgniter {-> handler.perform(true) }
        }
        else {
            handler.perform(true)
        }
        return channel
    }

}
