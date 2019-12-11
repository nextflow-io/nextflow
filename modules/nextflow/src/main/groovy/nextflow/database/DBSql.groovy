package nextflow.database

import groovy.util.logging.Slf4j
import groovy.transform.CompileStatic
import groovy.sql.Sql
import java.sql.ResultSetMetaData
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.extension.CH
import nextflow.Global
import nextflow.Session

@Slf4j
@CompileStatic
class DBSql {
    private DataflowWriteChannel target
    private String dbUrl
    private String dbUser
    private String dbPassword
    private String dbDriver
    private Sql connection
    private boolean delayConnection = false
    private String query
    private int columnCount = 0
    private List<String> columns = [] 
    private Closure metaclosure = { ResultSetMetaData meta ->
	columns = []
	columnCount = meta.getColumnCount()
        for (int i = 1; i <= columnCount; i++) {
	    columns << meta.getColumnLabel(i)
	}
    }

    protected void init(Map opts) {
        if( opts.dbUrl )
            dbUrl = opts.dbUrl
        if( opts.dbUser != null )
            dbUser = opts.dbUser
        if( opts.dbPassword )
            dbPassword = opts.dbPassword
        if (opts.dbDriver)
            dbDriver = opts.dbDriver
        if (opts.dbDelayConnecton != null)
            delayConnection = opts.dbDelayConnection as boolean
    }

    DBSql() {

    }

    DBSql(DataflowWriteChannel target, Map opts) {
        this.target = target
	if (opts.datasource) {
	    def session = Global.session as Session
	    def String configPath = "datasources.${opts.datasource}"
	    //def datasource = session.config.navigate("datasource.${opts.datasource}") //TODO better controls (build specific method
	    // TODO create a dedicated method to validate configuration
	    def Map datasource = [dbUrl: session.config.navigate("${configPath}.dbUrl") as String,
				  dbUser: session.config.navigate("${configPath}.dbUser") as String,
				  dbPassword: session.config.navigate("${configPath}.dbPassword") as String,
				  dbDriver: session.config.navigate("${configPath}.dbDriver") as String]
	    if (datasource){
		init(datasource)
	    } else {
		throw new IllegalArgumentException("Not a valid datasource definition: ${opts.datasource}")
	    }
	} else {
	    throw new IllegalArgumentException("Not enough information to establish a connection to the database")
	}
        if (!delayConnection)
            connect()
    }

    DBSql setQuery(String query) {
        this.query = query
        return this
    }

    void connect() {
        this.connection=Sql.newInstance(dbUrl, dbUser, dbPassword, dbDriver)
    }


    DataflowWriteChannel apply(Closure closure) {
        if( target == null )
            target = new DataflowQueue()
	    
        if( connection == null )
            connect()
	    
        query0(query, closure)
        
        target.bind(Channel.STOP)
        connection.close()
        return target
    }

    protected void query0(String query, Closure closure) {
        connection.eachRow(query, metaclosure) { row ->
	    Map<String, Object> lhm = new LinkedHashMap<String, Object>(columnCount, 1);
	    for (int i = 1; i <= columnCount; i++) {
		lhm.put(columns.get(i-1), row.getObject(i));
	    }
            target.bind(closure.call(lhm))
        }
    }
}

