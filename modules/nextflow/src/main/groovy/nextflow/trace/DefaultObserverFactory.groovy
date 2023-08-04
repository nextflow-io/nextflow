package nextflow.trace

import java.nio.file.Path

import nextflow.NextflowMeta
import nextflow.Session
import nextflow.script.ScriptBinding

/**
 * Creates Nextflow observes object
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultObserverFactory implements TraceObserverFactory {

    private Map config
    private Session session
    private ScriptBinding binding

    @Override
    Collection<TraceObserver> create(Session session) {
        this.session = session
        this.config = session.config
        this.binding = createBinding(session)

        final result = new ArrayList(10)
        createTraceFileObserver(result)
        createReportObserver(result)
        createTimelineObserver(result)
        createDagObserver(result)
        createAnsiLogObserver(result)
        return result
    }

    protected ScriptBinding createBinding(Session session) {
        final binding = new ScriptBinding( session.binding.getVariables() )
        binding.setVariable( 'workflow', session.workflowMetadata )
        binding.setVariable( 'nextflow', NextflowMeta.instance )

        return binding
    }

    protected void createAnsiLogObserver(Collection<TraceObserver> result) {
        if( session.ansiLog ) {
            session.ansiLogObserver = new AnsiLogObserver()
            result << session.ansiLogObserver
        }
    }

    /**
     * Create workflow report file observer
     */
    protected void createReportObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('report.enabled') as Boolean
        if( !isEnabled )
            return

        String fileName = config.navigateDynamic('report.file', binding)
        def maxTasks = config.navigate('report.maxTasks', ReportObserver.DEF_MAX_TASKS) as int
        if( !fileName ) fileName = ReportObserver.DEF_FILE_NAME
        def report = (fileName as Path).complete()
        def observer = new ReportObserver(report)
        observer.maxTasks = maxTasks
        config.navigate('report.overwrite') { observer.overwrite = it }
        result << observer
    }

    /**
     * Create timeline report file observer
     */
    protected void createTimelineObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('timeline.enabled') as Boolean
        if( !isEnabled )
            return

        String fileName = config.navigateDynamic('timeline.file', binding)
        if( !fileName ) fileName = TimelineObserver.DEF_FILE_NAME
        def traceFile = (fileName as Path).complete()
        def observer = new TimelineObserver(traceFile)
        config.navigate('timeline.overwrite')  { observer.overwrite = it }
        result << observer
    }

    protected void createDagObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('dag.enabled') as Boolean
        if( !isEnabled )
            return

        String fileName = config.navigateDynamic('dag.file', binding)
        if( !fileName ) fileName = GraphObserver.DEF_FILE_NAME
        def traceFile = (fileName as Path).complete()
        def observer = new GraphObserver(traceFile)
        config.navigate('dag.overwrite')  { observer.overwrite = it }
        result << observer
    }

    /*
     * create the execution trace observer
     */
    protected void createTraceFileObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('trace.enabled') as Boolean
        if( !isEnabled )
            return

        String fileName = config.navigateDynamic('trace.file', binding)
        if( !fileName ) fileName = TraceFileObserver.DEF_FILE_NAME
        def traceFile = (fileName as Path).complete()
        def observer = new TraceFileObserver(traceFile)
        config.navigate('trace.raw') { it -> observer.useRawNumbers(it == true) }
        config.navigate('trace.sep') { observer.separator = it }
        config.navigate('trace.fields') { observer.setFieldsAndFormats(it) }
        config.navigate('trace.overwrite') { observer.overwrite = it }
        result << observer
    }

}
