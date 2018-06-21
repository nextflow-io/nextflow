package nextflow.trace

import nextflow.Session
import spock.lang.Specification

/**
 * Created by edgar on 19/06/18.
 */
class ProvObserverTest extends Specification{

    def 'get Author Information from Manifest' (){
            //get Author name
            //get orcid
        given:
        def manifest = [author:'AuthorFoo']
        def session = Mock(Session)
        session.getConfig() >> [manifest: manifest]

        def observer = new ProvObserver()

        when:
        def result = observer.getAuthor(session)
        then:
        result == 'AuthorFoo'


    }
    def 'get Manifest' (){
        //get manifest
    }
    def 'get information for manifest file' (){
        //bundle/.ro/manifest.xml
        //author
        //orcid
        //createdBy

    }
    def 'get container technology information' (){
        //technology
        // container name
        //container sha256
    }
    def 'get information for metadata file' (){
        //technology
        // container name
        //container sha256
        //command line
        // uuid
        //Nf version
    }
    def 'check input entities' (){

    }
    def 'check output entities' (){

    }
    def 'check bundle manifest' (){

    }
}
