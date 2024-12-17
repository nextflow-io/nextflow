package nextflow.processor

import nextflow.Session
import nextflow.script.BaseScript
import nextflow.script.ScriptBinding
import nextflow.script.ScriptFile
import spock.lang.Specification

import test.TestHelper


class PublishOffloadManagerTest extends Specification {

    def 'should create task processor'(){
        given:
        def session = new Session();
        def scriptBinding = new ScriptBinding(session: session)
        def script = Stub(BaseScript)
        script.getBinding() >> scriptBinding
        def folder = TestHelper.createInMemTempDir()
        def file = folder.resolve('pipeline.nf'); file.text = 'println "hello"'
        def scriptFile = new ScriptFile(file)
        session.init(scriptFile)
        //session.start()
        session.script = script;
        def poManager = new PublishOffloadManager(session);
        when:
        poManager.init()
        then:
        poManager.copyProcessor != null
        poManager.moveProcessor != null
        cleanup:
        session.classesDir?.deleteDir()


    }

}
