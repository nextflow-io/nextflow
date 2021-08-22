package nextflow.sql

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.NF
import nextflow.Session
import nextflow.extension.CH
import nextflow.extension.ChannelExtensionPoint
import nextflow.extension.DataflowHelper
import nextflow.plugin.Scoped
import nextflow.sql.config.SqlConfig
import nextflow.sql.config.SqlDatasource
import nextflow.util.CheckHelper
/**
 * Provide a channel factory extension that allows the execution of Sql queries
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Scoped('sql')
class ChannelSqlExtension extends ChannelExtensionPoint {

    private static final Map QUERY_PARAMS = [dataSource: String]

    private static final Map INSERT_PARAMS = [dataSource: CharSequence, into: CharSequence, columns: [CharSequence, List], statement: CharSequence]

    private Session session
    private SqlConfig config

    protected void init(Session session) {
        this.session = session
        this.config = new SqlConfig((Map) session.config.dataSources)
    }

    DataflowWriteChannel fromQuery(String query) {
        fromQuery(Collections.emptyMap(), query)
    }

    DataflowWriteChannel fromQuery(Map opts, String query) {
        CheckHelper.checkParams('fromQuery', opts, QUERY_PARAMS)
        return queryToChannel(query, opts)
    }

    protected DataflowWriteChannel queryToChannel(String query, Map opts) {
        final channel = CH.create()
        final dataSource = dataSourceFromOpts(opts)
        final handler = new QueryHandler().withDatasource(dataSource).withStatement(query).withTarget(channel)
        if(NF.dsl2) {
            session.addIgniter {-> handler.perform(true) }
        }
        else {
            handler.perform(true)
        }
        return channel
    }

    protected SqlDatasource dataSourceFromOpts(Map opts) {
        final dsName = (opts?.dataSource ?: 'default') as String
        final dataSource = config.getDatasource(dsName)
        if( dataSource==null )
            throw new IllegalArgumentException("Unknown dataSource name: $dsName")
        return dataSource
    }

    @CompileDynamic
    DataflowWriteChannel sqlInsert( DataflowReadChannel source, Map opts=null ) {
        CheckHelper.checkParams('sqlInsert', opts, INSERT_PARAMS)
        final dataSource = dataSourceFromOpts(opts)
        final target = CH.createBy(source)
        final stopOnFirst = source instanceof DataflowExpression
        final insert = new InsertHandler(dataSource, opts)

        DataflowHelper.newOperator(source, target) { it ->
            // perform sql insert
            insert.perform(it)
            // bind to the target
            target.bind(it)
            // check for termination
            if( it == Channel.STOP || stopOnFirst ) {
                ((DataflowProcessor) getDelegate()).terminate()
                insert.close()
            }
        }

        return target
    }

}
