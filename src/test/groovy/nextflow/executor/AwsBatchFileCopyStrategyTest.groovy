package nextflow.executor
import java.nio.file.Paths

import nextflow.processor.TaskBean
import spock.lang.Specification
import test.TestHelper

class AwsBatchFileCopyStrategyTest extends Specification {

    def 'should strip out file/folder name from target S3 path' () {
        given:
        def bean = new TaskBean(
          inputFiles: [one:'bam1',two:'bam2'],
          workDir: Paths.get('/data/work/a4/a4cf1fd2b040ae9b36ff46c34'),
          targetDir: Paths.get('/data/results'),
          outputFiles: ["outputs_*","final_folder"]
        )
        def fileCopy = new AwsBatchFileCopyStrategy(bean, Mock(AwsOptions))
        expect:
        fileCopy.touchFile(Paths.get("/data/work/41/41393a4cf1fd2b040ae9b36ff46c34/.command.run")) == "touch .command.run && nxf_s3_upload .command.run s3://data/work/41/41393a4cf1fd2b040ae9b36ff46c34 && rm .command.run"
        fileCopy.copyFile("nobel_prize_results.gz",Paths.get("/data/work/41/41393a4cf1fd2b040ae9b36ff46c34/nobel_prize_results.gz")) == "nxf_s3_upload nobel_prize_results.gz s3://data/work/41/41393a4cf1fd2b040ae9b36ff46c34"
        fileCopy.exitFile(Paths.get("/data/work/41/41393a4cf1fd2b040ae9b36ff46c34/.exitcode")) == ".exitcode && nxf_s3_upload .exitcode s3://data/work/41/41393a4cf1fd2b040ae9b36ff46c34 || true"
        fileCopy.getUnstageOutputFilesScript() == "\nnxf_s3_upload 'outputs_*' s3://data/results || true\nnxf_s3_upload 'final_folder' s3://data/results || true"
    }

    def 'should check the beforeScript' () {

        given:
        def bean = new TaskBean(
                inputFiles: [:],
                workDir: Paths.get('/data/work/a4/a4cf1fd2b040ae9b36ff46c34'),
                outputFiles: []
        )
        def opts = Spy(AwsOptions)
        def copy = new AwsBatchFileCopyStrategy(bean, opts)

        when:
        def script = copy.getBeforeStartScript()
        then:
        1 * opts.getCliPath() >> null
        1 * opts.getStorageClass() >> null
        1 * opts.getStorageEncryption() >> null

        script == '''
                # aws helper
                nxf_s3_upload() {
                    local pattern=$1
                    local s3path=$2
                    for name in $pattern;do
                      if [[ -d "$name" ]]; then
                        aws s3 cp $name $s3path/$name --quiet --recursive --storage-class STANDARD
                      else
                        aws s3 cp $name $s3path/$name --quiet --storage-class STANDARD
                      fi
                  done
                }

            '''.stripIndent()

        when:
        script = copy.getBeforeStartScript()
        then:
        1 * opts.getCliPath() >> '/foo/aws'
        1 * opts.getStorageClass() >> 'REDUCED_REDUNDANCY'
        2 * opts.getStorageEncryption() >> 'AES256'

        script == '''
                # aws helper
                nxf_s3_upload() {
                    local pattern=$1
                    local s3path=$2
                    for name in $pattern;do
                      if [[ -d "$name" ]]; then
                        /foo/aws s3 cp $name $s3path/$name --quiet --recursive --sse AES256 --storage-class REDUCED_REDUNDANCY
                      else
                        /foo/aws s3 cp $name $s3path/$name --quiet --sse AES256 --storage-class REDUCED_REDUNDANCY
                      fi
                  done
                }

            '''.stripIndent()

    }


    def 'should return stage input input file'() {
        given:
        def file = TestHelper.createInMemTempFile('foo.txt')
        def folder = TestHelper.createInMemTempDir()

        def bean = new TaskBean(
                inputFiles: [:],
                workDir: Paths.get('/data/work/a4/a4cf1fd2b040ae9b36ff46c34'),
                outputFiles: []
        )
        def opts = Spy(AwsOptions)
        def copy = new AwsBatchFileCopyStrategy(bean, opts)

        when:
        def script = copy.stageInputFile( file, 'bar.txt')
        then:
        1 * opts.getCliPath() >> null
        script == "aws s3 cp --quiet s3:/$file bar.txt" as String

        when:
        script = copy.stageInputFile( folder, 'bar')
        then:
        1 * opts.getCliPath() >> null
        script == "aws s3 cp --quiet --recursive s3:/$folder bar" as String


        when:
        script = copy.stageInputFile( folder, 'bar')
        then:
        1 * opts.getCliPath() >> '/home/bin/aws'
        script == "/home/bin/aws s3 cp --quiet --recursive s3:/$folder bar" as String
    }
    
}
