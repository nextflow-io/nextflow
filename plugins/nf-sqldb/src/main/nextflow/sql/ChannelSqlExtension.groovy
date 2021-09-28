package nextflow.sql


import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import nextflow.Channel
import nextflow.NF
import nextflow.Session
import nextflow.extension.CH
import nextflow.extension.ChannelExtensionPoint
import nextflow.extension.DataflowHelper
import nextflow.plugin.Scoped
import nextflow.sql.config.SqlConfig
import nextflow.sql.config.SqlDataSource
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

    private static final Map QUERY_PARAMS = [db: CharSequence]

    private static final Map INSERT_PARAMS = [db: CharSequence, into: CharSequence, columns: [CharSequence, List], statement: CharSequence]

    private Session session
    private SqlConfig config

    protected void init(Session session) {
        this.session = session
        this.config = new SqlConfig((Map) session.config.navigate('sql.db'))
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
        final handler = new QueryHandler().withDataSource(dataSource).withStatement(query).withTarget(channel)
        if(NF.dsl2) {
            session.addIgniter {-> handler.perform(true) }
        }
        else {
            handler.perform(true)
        }
        return channel
    }

    protected SqlDataSource dataSourceFromOpts(Map opts) {
        final dsName = (opts?.db ?: 'default') as String
        final dataSource = config.getDataSource(dsName)
        if( dataSource==null ) {
            def msg = "Unknown db name: $dsName"
            def choices = config.getDataSourceNames().closest(dsName) ?: config.getDataSourceNames()
            if( choices?.size() == 1 )
                msg += " - Did you mean: ${choices.get(0)}?"
            else if( choices )
                msg += " - Did you mean any of these?\n" + choices.collect { "  $it"}.join('\n') + '\n'
            throw new IllegalArgumentException(msg)
        }
        return dataSource
    }

    DataflowWriteChannel sqlInsert( DataflowReadChannel source, Map opts=null ) {
        CheckHelper.checkParams('sqlInsert', opts, INSERT_PARAMS)
        final dataSource = dataSourceFromOpts(opts)
        final target = CH.createBy(source)
        final singleton = target instanceof DataflowExpression
        final insert = new InsertHandler(dataSource, opts)

        final next = { it ->
            insert.perform(it)
            target.bind(it)
        }

        final done = {
            insert.close()
            if( !singleton ) target.bind(Channel.STOP)
        }

        DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
        return target
    }

}
