package nextflow.lineage.fs

import java.nio.channels.SeekableByteChannel
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileTime
import java.time.Instant

class LinIntermediatePath extends LinMetadataPath{

    LinIntermediatePath(LinFileSystem fs, String path) {
        super("", FileTime.from(Instant.now()), fs, path, null)
    }

    @Override
    InputStream newInputStream() {
        throw new UnsupportedOperationException()
    }
    @Override
    SeekableByteChannel newSeekableByteChannel(){
        throw new UnsupportedOperationException()
    }
    @Override
    <A extends BasicFileAttributes> A readAttributes(Class<A> type){
        return (A) new BasicFileAttributes() {
            @Override
            long size() { return 0 }

            @Override
            FileTime lastModifiedTime() { FileTime.from(Instant.now()) }

            @Override
            FileTime lastAccessTime() { FileTime.from(Instant.now()) }

            @Override
            FileTime creationTime() { FileTime.from(Instant.now()) }

            @Override
            boolean isRegularFile() { return false }

            @Override
            boolean isDirectory() { return true }

            @Override
            boolean isSymbolicLink() { return false }

            @Override
            boolean isOther() { return false }

            @Override
            Object fileKey() { return null }
        }
    }
}
