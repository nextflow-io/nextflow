package nextflow.sraql


import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.NF
import nextflow.Session
import nextflow.extension.CH
import nextflow.extension.ChannelExtensionPoint
import nextflow.plugin.Scoped
import nextflow.sraql.config.SraqlConfig
import nextflow.sraql.config.SraqlDataSource
import nextflow.util.CheckHelper
/**
 * Provide a channel factory extension that allows the execution of SRAQL queries
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
@Slf4j
@CompileStatic
@Scoped('sraql')
class ChannelSraqlExtension extends ChannelExtensionPoint {

    private static final Map QUERY_PARAMS = [source: CharSequence]
    private Session session
    private SraqlConfig config

    protected void init(Session session) {
        this.session = session
        this.config = new SraqlConfig((String) session.config.navigate('sraql.source'))
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
        final handler = new QueryHandler()
                .withDataSource(dataSource)
                .withStatement(query)
                .withTarget(channel)

        if(NF.dsl2) {
            session.addIgniter {-> handler.perform(true) }
        }
        else {
            handler.perform(true)
        }
        return channel
    }

    protected SraqlDataSource dataSourceFromOpts(Map opts) {
        final dsName = opts.source  as String
        final dataSource = config.getDataSource(dsName)
        return dataSource
    }
}
