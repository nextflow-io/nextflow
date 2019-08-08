package nextflow.trace

import java.nio.file.Path

import nextflow.Session

/**
 * Creates Nextflow observes object
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultObserverFactory implements TraceObserverFactory {

    private Map config
    private Session session

    @Override
    Collection<TraceObserver> create(Session session) {
        this.session = session
        this.config = session.config

        final result = new ArrayList(10)
        createTraceFileObserver(result)
        createReportObserver(result)
        createTimelineObserver(result)
        createDagObserver(result)
        createWebLogObserver(result)
        createAnsiLogObserver(result)
        return result
    }

    protected void createAnsiLogObserver(Collection<TraceObserver> result) {
        if( session.ansiLog ) {
            session.ansiLogObserver = new AnsiLogObserver()
            result << session.ansiLogObserver
        }
    }

    /**
     * Create workflow message observer
     * @param result
     */
    protected void createWebLogObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('weblog.enabled') as Boolean
        String url = config.navigate('weblog.url') as String
        if ( isEnabled ) {
            if ( !url ) url = WebLogObserver.DEF_URL
            def observer = new WebLogObserver(url)
            result << observer
        }
    }

    /**
     * Create workflow report file observer
     */
    protected void createReportObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('report.enabled') as Boolean
        if( !isEnabled )
            return

        String fileName = config.navigate('report.file')
        def maxTasks = config.navigate('report.maxTasks', ReportObserver.DEF_MAX_TASKS) as int
        if( !fileName ) fileName = ReportObserver.DEF_FILE_NAME
        def report = (fileName as Path).complete()
        def observer = new ReportObserver(report)
        observer.maxTasks = maxTasks

        result << observer
    }

    /**
     * Create timeline report file observer
     */
    protected void createTimelineObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('timeline.enabled') as Boolean
        if( !isEnabled )
            return

        String fileName = config.navigate('timeline.file')
        if( !fileName ) fileName = TimelineObserver.DEF_FILE_NAME
        def traceFile = (fileName as Path).complete()
        def observer = new TimelineObserver(traceFile)
        result << observer
    }

    protected void createDagObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('dag.enabled') as Boolean
        if( !isEnabled )
            return

        String fileName = config.navigate('dag.file')
        if( !fileName ) fileName = GraphObserver.DEF_FILE_NAME
        def traceFile = (fileName as Path).complete()
        def observer = new GraphObserver(traceFile)
        result << observer
    }

    /*
     * create the execution trace observer
     */
    protected void createTraceFileObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('trace.enabled') as Boolean
        if( !isEnabled )
            return

        String fileName = config.navigate('trace.file')
        if( !fileName ) fileName = TraceFileObserver.DEF_FILE_NAME
        def traceFile = (fileName as Path).complete()
        def observer = new TraceFileObserver(traceFile)
        config.navigate('trace.raw') { it -> observer.useRawNumbers(it == true) }
        config.navigate('trace.sep') { observer.separator = it }
        config.navigate('trace.fields') { observer.setFieldsAndFormats(it) }
        result << observer
    }

}
