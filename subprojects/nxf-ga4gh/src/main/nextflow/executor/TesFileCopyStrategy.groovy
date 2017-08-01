package nextflow.executor

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.util.Escape

@CompileStatic
class TesFileCopyStrategy implements ScriptFileCopyStrategy {

    @Override
    String getBeforeStartScript() {
        return ''
    }

    @Override
    String getStageInputFilesScript() {
        return ''
    }

    @Override
    String getUnstageOutputFilesScript() {
        return ''
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

    @Override
    String exitFile(Path file) {
        Escape.path(file.getFileName())
    }
}