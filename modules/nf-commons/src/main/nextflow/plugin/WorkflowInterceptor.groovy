package nextflow.plugin

import org.pf4j.ExtensionPoint

/**
 * Extension point for intercepting named workflow execution.
 *
 * Plugins implement this interface to gain control before
 * each named workflow runs. The interceptor uses the
 * {@code proceed} callback to decide whether to execute
 * the original workflow logic:
 * - Call {@code proceed()} to execute normally, optionally
 *   post-processing the result (e.g. archiving)
 * - Return without calling {@code proceed()} to skip
 *   execution (e.g. restore from archive)
 *
 * Entry workflows (name is null) are never intercepted.
 *
 * @author Zhibo Huang
 */
interface WorkflowInterceptor extends ExtensionPoint {

    /**
     * Intercept the execution of a named workflow.
     *
     * @param workflow the named workflow definition (WorkflowDef)
     * @param args the arguments passed to the workflow
     * @param proceed callback that executes the original workflow logic
     * @return the workflow execution result (ChannelOut)
     */
    Object intercept(Object workflow, Object[] args, Closure proceed)
}
