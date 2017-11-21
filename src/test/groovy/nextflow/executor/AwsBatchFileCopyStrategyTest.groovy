package nextflow.executor
import java.nio.file.Paths

import nextflow.processor.TaskBean
import spock.lang.Specification
import test.TestHelper

class AwsBatchFileCopyStrategyTest extends Specification {

    def 'should strip out file/folder name from target S3 path' () {
        given:
        def OUTPUTS = ["outputs_*","final_folder"]
        def TARGET =  Paths.get('/data/results')
        def FILE = Paths.get('/some/data/nobel_prize_results.gz')
        def EXIT = Paths.get('/some/path/.exitcode')
        def RUN = Paths.get('/some/data/.command.run')
        def copy = new AwsBatchFileCopyStrategy(Mock(TaskBean), Mock(AwsOptions))
        expect:
        copy.touchFile(RUN) == "echo start | aws s3 cp - s3://some/data/.command.run"
        copy.copyFile("nobel_prize_results.gz",Paths.get("/some/data/nobel_prize_results.gz")) == "nxf_s3_upload nobel_prize_results.gz s3://some/data"
        copy.exitFile(EXIT) == "| aws s3 cp - s3://some/path/.exitcode || true"
        copy.getUnstageOutputFilesScript(OUTPUTS,TARGET) == "\nnxf_s3_upload 'outputs_*' s3://data/results || true\nnxf_s3_upload 'final_folder' s3://data/results || true"
        copy.stageInputFile(FILE, 'foo.txt') == 'aws s3 cp --quiet s3://some/data/nobel_prize_results.gz foo.txt'
    }

    def 'should check the beforeScript' () {

        given:
        def bean = Mock(TaskBean)
        def opts = Mock(AwsOptions)
        def copy = Spy(AwsBatchFileCopyStrategy, constructorArgs: [bean, opts])

        when:
        def script = copy.getBeforeStartScript()
        then:
        1 * opts.getCliPath() >> null
        1 * opts.getStorageClass() >> null
        1 * opts.getStorageEncryption() >> null

        script ==   '''
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

    def 'should return env variables' () {

        given:
        def ENV = [FOO: 'hola', BAR:'world', PATH:'xxx']
        def bean = Mock(TaskBean)
        def opts = Mock(AwsOptions)
        def copy = Spy(AwsBatchFileCopyStrategy, constructorArgs: [bean, opts])

        when:
        def script = copy.getEnvScript(ENV)
        then:
        // note: PATH is always removed
        opts.getRemoteBinDir() >> null
        script == '''
            export BAR="world"
            export FOO="hola"

            '''.stripIndent().leftTrim()

        when:
        script = copy.getEnvScript(ENV)
        then:
        opts.getRemoteBinDir() >> '/foo/bar'
        script == '''
            aws s3 cp --recursive --quiet s3://foo/bar $PWD/nextflow-bin
            chmod +x $PWD/nextflow-bin/*
            export PATH=$PWD/nextflow-bin:$PATH
            export BAR="world"
            export FOO="hola"

            '''.stripIndent().leftTrim()

        when:
        script = copy.getEnvScript(ENV)
        then:
        opts.getCliPath() >> '/conda/bin/aws'
        opts.getRemoteBinDir() >> '/foo/bar'
        script == '''
            /conda/bin/aws s3 cp --recursive --quiet s3://foo/bar $PWD/nextflow-bin
            chmod +x $PWD/nextflow-bin/*
            export PATH=$PWD/nextflow-bin:$PATH
            export BAR="world"
            export FOO="hola"

            '''.stripIndent().leftTrim()

        when:
        script = copy.getEnvScript(ENV)
        then:
        opts.getCliPath() >> '/conda/bin/aws'
        opts.getRemoteBinDir() >> '/foo/bar'
        opts.getRegion() >> 'eu-west-1'
        script == '''
            /conda/bin/aws s3 cp --recursive --quiet s3://foo/bar $PWD/nextflow-bin
            chmod +x $PWD/nextflow-bin/*
            export PATH=$PWD/nextflow-bin:$PATH
            export BAR="world"
            export FOO="hola"

            '''.stripIndent().leftTrim()

    }


    def 'should return stage input input file'() {
        given:
        def file = TestHelper.createInMemTempFile('foo.txt')
        def folder = TestHelper.createInMemTempDir()

        def bean = Mock(TaskBean)
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
