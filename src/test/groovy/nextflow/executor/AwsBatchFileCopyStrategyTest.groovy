package nextflow.executor

import spock.lang.Specification
import java.nio.file.*
import nextflow.processor.TaskBean

class AwsBatchFileCopyStrategyTest extends Specification {

    def 'should strip out file/folder name from target S3 path' () {
        given:
        def bean = new TaskBean(
          inputFiles: [one:'bam1',two:'bam2'],
          workDir: Paths.get('/data/work/a4/a4cf1fd2b040ae9b36ff46c34'),
          targetDir: Paths.get('/data/results'),
          outputFiles: ["outputs_*","final_folder"]
        )
        def fileCopy = new AwsBatchFileCopyStrategy(bean)
        expect:
        fileCopy.touchFile(Paths.get("/data/work/41/41393a4cf1fd2b040ae9b36ff46c34/.command.run")) == "touch .command.run && nxf_s3_upload .command.run s3://data/work/41/41393a4cf1fd2b040ae9b36ff46c34 && rm .command.run"
        fileCopy.copyFile("nobel_prize_results.gz",Paths.get("/data/work/41/41393a4cf1fd2b040ae9b36ff46c34/nobel_prize_results.gz")) == "nxf_s3_upload nobel_prize_results.gz s3://data/work/41/41393a4cf1fd2b040ae9b36ff46c34"
        fileCopy.exitFile(Paths.get("/data/work/41/41393a4cf1fd2b040ae9b36ff46c34/.exitcode")) == ".exitcode && nxf_s3_upload .exitcode s3://data/work/41/41393a4cf1fd2b040ae9b36ff46c34 || true"
        fileCopy.getUnstageOutputFilesScript() == "\nnxf_s3_upload 'outputs_*' s3://data/results || true\nnxf_s3_upload 'final_folder' s3://data/results || true"
    }
    
}
