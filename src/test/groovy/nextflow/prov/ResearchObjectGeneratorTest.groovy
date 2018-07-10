package nextflow.prov

import nextflow.Session
import nextflow.trace.ProvObserver
import spock.lang.Specification

/**
 * Created by edgar on 21/06/18.
 */
//      how to declare a mock
//def session = Mock(Session)
//session.getConfig() >> [manifest: manifest]
class ResearchObjectGeneratorTest extends Specification {

    def 'check author from Manifest' (){
        given:
        def manifest = [author:'AuthorFoo']
        def observer = new ResearchObjectGenerator()
        def bundle = observer.generateROBundle()

        when:
        observer.setManifest(bundle,manifest)
        def resultManifest = bundle.getManifest()

        then:
        resultManifest.getCreatedBy().getName() == 'AuthorFoo'
    }
    def 'check null Author from Manifest' (){
        given:
        def manifest = [orcid:'ORCIDFoo']
        def observer = new ResearchObjectGenerator()
        def bundle = observer.generateROBundle()

        when:
        observer.setManifest(bundle,manifest)
        def resultManifest = bundle.getManifest()

        then:
        resultManifest.getCreatedBy()== null
    }
    def 'check ORCID from Manifest' (){
        given:
        def manifest = [author:'AuthorFoo', ORCID:'ORCIDFoo']
        def observer = new ResearchObjectGenerator()
        def bundle = observer.generateROBundle()

        when:
        observer.setManifest(bundle,manifest)
        def resultManifest = bundle.getManifest()

        then:
        resultManifest.getCreatedBy().getOrcid().toString() == 'ORCIDFoo'
    }
    def 'check null ORCID from Manifest' (){
        given:
        def manifest = [author:'AuthorFoo']
        def observer = new ResearchObjectGenerator()
        def bundle = observer.generateROBundle()

        when:
        observer.setManifest(bundle,manifest)
        def resultManifest = bundle.getManifest()

        then:
        resultManifest.getCreatedBy().getOrcid().toString() == '**ORCID_not_provided**'
    }
    def 'check null author name, means NOT ORCID'(){
        given:
        def manifest = [ ORCID:'ORCIDFoo']
        def observer = new ResearchObjectGenerator()
        def bundle = observer.generateROBundle()

        when:
        observer.setManifest(bundle,manifest)
        def resultManifest = bundle.getManifest()

        then:
        resultManifest.getCreatedBy() == null
    }
    def 'check Data folder'(){
        //check all input data is inside
        //check what happens when no input is defined
    }

    def 'check Workflow folder'(){
        //check all baseDir files are inside
    }

    def 'check Snapshot folder'(){
        //check all baseDir files are inside
        //check is the .command.sh is there
    }
    def 'check Metadata folder'(){
        //check metadata.json file

    }
    def 'check log file'(){
        //check the correctness of the info
        //check what happen when NO author
        //check what happen wih NO container
    }

    def 'check provenance file'(){
        //check if prov file exist
    }
}
