package nextflow.cloud.azure.batch

import java.nio.file.Path

import nextflow.executor.SimpleFileCopyStrategy
import nextflow.processor.TaskBean

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzFileCopyStrategy extends SimpleFileCopyStrategy {

    protected AzFileCopyStrategy() {}

    AzFileCopyStrategy(TaskBean bean) {
        super(bean)
    }

    @Override
    String getStageInputFilesScript(Map<String, Path> inputFiles) {
        // null
        return null
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String stageInputFile( Path path, String targetName ) {
        // todo
        return null
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {
        // todo
        null
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String touchFile( Path file ) {
        super.touchFile( workDir.relativize(file) )
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String fileStr( Path path ) {
        super.fileStr( workDir.relativize(path) )
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String copyFile( String name, Path target ) {
        super.copyFile(name, workDir.relativize(target))
    }

    /**
     * {@inheritDoc}
     */
    String exitFile( Path path ) {
        super.exitFile( workDir.relativize(path) )
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String pipeInputFile( Path path ) {
        super.pipeInputFile( workDir.relativize(path) )
    }

}
