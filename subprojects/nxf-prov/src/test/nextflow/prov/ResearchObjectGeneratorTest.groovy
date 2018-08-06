package nextflow.prov

import org.apache.taverna.robundle.Bundle
import spock.lang.Specification
import java.nio.file.Path
import java.nio.file.Paths

/**
 * Created by edgar on 21/06/18.
 */
class ResearchObjectGeneratorTest extends Specification {

    def 'check generateFileStructure()'(){
        given:
        def observer = Spy(ResearchObjectGenerator)
        def bundle = Mock(Bundle)
        def metaFolderPath = Mock(Path)
        def actualPath = Paths.get(System.getProperty("user.dir"))

        when:
        bundle.getRoot() >> actualPath
        observer.generateFileStructure(bundle)

        then:
        //1 * bundle.getRoot().resolve("metadata") >> metaFolderPath
        1 * observer.createROdirectory(Paths.get("${actualPath}/outputs"))>> null
        1 * observer.createROdirectory(Paths.get("${actualPath}/data"))>> null
        1 * observer.createROdirectory(Paths.get("${actualPath}/workflow"))>> null
        1 * observer.createROdirectory(Paths.get("${actualPath}/metadata"))>> null
        1 * observer.createROdirectory(Paths.get("${actualPath}/snapshot"))>> null
    }

    def 'check setManifest()' (){
        given:
        def observer = Spy(ResearchObjectGenerator)
        def bundle = Mock(Bundle)
        def manifest = Mock(org.apache.taverna.robundle.manifest.Manifest)

        when:
        observer.setManifest(bundle)

        then:
        1 * bundle.getManifest() >> manifest
        1 * manifest.setId(URI.create('/')) >> null
        1 * observer.setAuthorInformation(manifest) >> null
        1 * observer.setAggregationManifest(manifest) >> null
    }

    def 'check generateLogFile()'(){
        given:
        def observer = Spy(ResearchObjectGenerator)
        def bundle = Mock(Bundle)
        def logFile = Mock(File)
        def actualPath = Paths.get(System.getProperty("user.dir"))
        logFile.name >>"log.txt"
        logFile.path >> actualPath

        when:
        bundle.getRoot() >> actualPath
        observer.generateFileStructure(bundle)
        observer.generateLogFile(bundle)

        then:
        1 * observer.getLogInfo() >> logFile
        1 * observer.fileToBundle(bundle,actualPath,logFile.name, "metadata")
    }

   /* def 'check addProvenanceFile()'(){
        given:
        def observer = Spy(ResearchObjectGenerator)
        def bundle = Mock(Bundle)
        def provFile = Mock(File)
        def actualPath = Paths.get(System.getProperty("user.dir"))
        provFile.name >>"provenance.json"
        provFile.path >> actualPath

        when:
        bundle.getRoot() >> actualPath
        observer.generateFileStructure(bundle)
        observer.addProvenanceFile(bundle)

        then:
        1 * observer.fileToBundle(bundle,actualPath,provFile.name, "metadata")
    }*/
/*
    def 'check generateDataFolder()'(){
        given:
        def observer = Spy(ResearchObjectGenerator)
        def bundle = Mock(Bundle)

        when:
        def numberFiles = observer.inputFiles.unique().size()

        then:
        numberFiles * observer.fileToBundle(bundle, auxPath, auxPath.getFileName().toString(), dataFolderName)
    }*/


    def 'check null Author into RO bundle' (){
        given:
        def observer = new ResearchObjectGenerator()
        def bundle = observer.generateROBundle()

        when:
        observer.setManifest(bundle)
        def resultManifest = bundle.getManifest()

        then:
        resultManifest.getCreatedBy()== null
    }

    def 'check ORCID into RO bundle' (){
        given:
        def observer = new ResearchObjectGenerator()
        observer.author='AuthorFoo'
        observer.authorORCID='ORCIDFoo'
        def bundle = observer.generateROBundle()

        when:
        observer.setManifest(bundle)
        def resultManifest = bundle.getManifest()

        then:
        resultManifest.getCreatedBy().getOrcid().toString() == 'ORCIDFoo'
    }
    def 'check null ORCID ' (){
        given:
        def observer = new ResearchObjectGenerator()
        observer.author='AuthorFoo'
        def bundle = observer.generateROBundle()

        when:
        observer.setManifest(bundle)
        def resultManifest = bundle.getManifest()

        then:
        resultManifest.getCreatedBy().getOrcid().toString() == '**ORCID_not_provided**'
    }

    def 'check null author name, means NOT ORCID'(){
        given:
        def observer = new ResearchObjectGenerator()
        observer.authorORCID='ORCIDFoo'
        def bundle = observer.generateROBundle()

        when:
        observer.setManifest(bundle)
        def resultManifest = bundle.getManifest()

        then:
        resultManifest.getCreatedBy() == null
    }
    /*
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
 */
}
