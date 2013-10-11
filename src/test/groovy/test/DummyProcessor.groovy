package test

import nextflow.Session
import nextflow.executor.NopeExecutor
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.script.BaseScript

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DummyProcessor extends TaskProcessor {

    DummyProcessor(Session session, BaseScript script, TaskConfig taskConfig) {
        super(new NopeExecutor(), session, new DummyScript(), taskConfig, {})
    }

    @Override
    protected void createOperator() {
    }
}



class DummyScript extends BaseScript {
    @Override
    Object run() {
        return null
    }

}