package nextflow.prov

import org.apache.taverna.robundle.Bundle
import spock.lang.Specification
import java.nio.file.Path
import java.nio.file.Paths

/**
 *
 * @author Edgar Garriga <edgano@gmail.com>
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
        //1 * bundle.getRoot().resolve("metadata") >> Path
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
        1 * observer.setLogInfo() >> logFile
        1 * observer.fileToBundle(bundle,actualPath,logFile.name, "metadata")
    }

    def 'check addProvenanceFile()'(){
        given:
        def observer = Spy(ResearchObjectGenerator)
        def bundle = Mock(Bundle)
        def provFile = Mock(File)
        def actualPath = Paths.get(System.getProperty("user.dir"))
        provFile.name >>"provenance.json"
        provFile.path >>"provenance.json"

        when:
        bundle.getRoot() >> actualPath
        observer.generateFileStructure(bundle)
        observer.addProvenanceFile(bundle)

        then:
        1 * observer.fileToBundle(bundle,Paths.get(provFile.path),provFile.name, "metadata") >> null
    }

    def 'check generateDataFolder()'(){
        given:
        def observer = Spy(ResearchObjectGenerator)
        def bundle = Mock(Bundle)
        def provFile = Mock(File)
        def actualPath = Paths.get(System.getProperty("user.dir"))
        def numberFiles = observer.inputFiles.unique().size()
        def dataFolderName = "data"

        when:
        bundle.getRoot() >> actualPath
        observer.generateDataFolder()

        then:
        numberFiles * observer.fileToBundle(bundle, provFile.getPath(), provFile.getName(), dataFolderName)
    }

    def 'check generateOutputFolder()'(){
        given:
        def observer = Spy(ResearchObjectGenerator)
        def bundle = Mock(Bundle)
        def provFile = Mock(File)
        def actualPath = Paths.get(System.getProperty("user.dir"))
        def numberFiles = observer.outputFiles.unique().size()
        def outputFolderName = "outputs"

        when:
        bundle.getRoot() >> actualPath
        observer.generateOutputFolder()

        then:
        numberFiles * observer.fileToBundle(bundle, provFile.getPath(), provFile.getName(), outputFolderName)    }


    def 'check generateWorkflowFolder()'(){
        given:
        def observer = Spy(ResearchObjectGenerator)
        def bundle = Mock(Bundle)
        def actualPath = Paths.get(System.getProperty("user.dir"))
        observer.baseDir = "${actualPath}/test"

        when:
        bundle.getRoot() >> actualPath
        observer.generateWorkflowFolder(bundle)

        then:
        1 * observer.getFilesFromDir(observer.baseDir)
    }

    def 'check generateSnapshot'(){
        given:
        def observer = Spy(ResearchObjectGenerator)
        def bundle = Mock(Bundle)
        def actualPath = Paths.get(System.getProperty("user.dir"))
        observer.commandLine = "holaMundo"

        when:
        bundle.getRoot() >> actualPath
        observer.generateSnapshot(bundle)

        then:
        1 * observer.generateScript()
    }

    def 'check generateMetadataFolder'(){
        given:
        def observer = Spy(ResearchObjectGenerator)
        def bundle = Mock(Bundle)
        Path actualPath = Paths.get(System.getProperty("user.dir"))

        when:
        bundle.getRoot() >> actualPath
        observer.generateMetadataFolder(bundle)

        then:
        1 * observer.generateMetadataFile()
        //1 * observer.fileToBundle(bundle, metadataFile,metadataFileName,metadataFolder)

    }

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
}
