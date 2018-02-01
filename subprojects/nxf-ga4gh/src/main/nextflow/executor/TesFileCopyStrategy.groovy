package nextflow.executor

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.util.Escape

@CompileStatic
class TesFileCopyStrategy implements ScriptFileCopyStrategy {

    @Override
    String getBeforeStartScript() {
        return null
    }

    @Override
    Map<String, Path> resolveForeignFiles(Map<String, Path> inputFile) {
        return inputFile
    }

    @Override
    String getStageInputFilesScript(Map<String, Path> inputFiles) {
        return null
    }

    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {
        return null
    }

    @Override
    String touchFile(Path file) {
        return ''
    }

    @Override
    String fileStr(Path file) {
        Escape.path(file.getFileName())
    }

    @Override
    String copyFile(String name, Path target) {
        return 'true'
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String exitFile(Path file) {
        "> ${Escape.path(file.getName())}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String pipeInputFile(Path file) {
        " < ${Escape.path(file.getName())}"
    }

    @Override
    String getEnvScript(Map environment, String wrapName) {
        // TODO this must be implement as SimpleFileCopyStrategy#getEnvScript
        return null
    }
}