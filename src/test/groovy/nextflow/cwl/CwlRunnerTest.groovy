package nextflow.cwl
import java.nio.file.Paths

import nextflow.script.BaseScript
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.ValueInParam
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CwlRunnerTest extends Specification {

    def 'should parse cwl command line tool' () {
        given:
        def text = '''
        cwlVersion: v1.0
        class: CommandLineTool
        baseCommand: samtools
        arguments: [index]

        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - $(inputs.bamfile)

        hints:
          - class: ResourceRequirement
            coresMin: 2
            outdirMin: 2048
            ramMin: 4096

        inputs:
          bamfile:
            type: File
            inputBinding:
              position: 1
          outbam:
            type: string
            default: test.bai
            inputBinding:
              position: 2

        outputs:
          indexout:
            type: File[]
            outputBinding:
              glob: "*.bam"
        '''
                .stripIndent()

        def INPUTS = [
                bamfile: Paths.get('/some/path'),
                outbam: 'foo'
        ]

        def cwl = new CwlRunner(Mock(BaseScript), INPUTS)

        when:
        cwl.parse(text)

        then:
        // the resulting command
        cwl.command == 'samtools index ${bamfile} ${outbam}'
        // process inputs
        cwl.config.getInputs().get(0) instanceof FileInParam
        cwl.config.getInputs().get(0).getName() == 'bamfile'
        cwl.config.getInputs().get(1) instanceof ValueInParam
        cwl.config.getInputs().get(1).getName() == 'outbam'
        // process outputs
        cwl.config.getOutputs().get(0) instanceof FileOutParam
        cwl.config.getOutputs().get(0).getFilePattern() == '*.bam'

    }
    def 'parse baseCommand list' () {
        given:
        def text = '''
        cwlVersion: v1.0
        class: CommandLineTool
        baseCommand: tar
        inputs:
          tarfile:
            type: File
            inputBinding:
              position: 1
              prefix: xf
          extractfile:
            type: string
            inputBinding:
              position: 2
        outputs:
          example_out:
            type: File
            outputBinding:
              glob: $(inputs.extractfile)
        '''
                .stripIndent()

        def INPUTS = [
                tarfile: Paths.get('/some/path'),
                extractfile: 'Hello.java'
        ]

        def cwl = new CwlRunner(Mock(BaseScript), INPUTS)

        when:
        cwl.parse(text)

        then:
        // the resulting command
        cwl.command == 'tar xf ${tarfile} ${extractfile}'

    }
}
