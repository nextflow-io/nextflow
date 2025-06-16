package software.amazon.nio.spi.s3

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.TagAwareFile
import nextflow.util.TestOnly
import software.amazon.awssdk.services.s3.model.Tag

import java.nio.file.FileSystem
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.WatchEvent
import java.nio.file.WatchKey
import java.nio.file.WatchService

@Slf4j
@CompileStatic
class NextflowS3Path implements Path, TagAwareFile{

    // Delegate annotation was not working due to package private
    private final S3Path delegate

    private Map<String,String> tags

    private String contentType

    private String storageClass

    NextflowS3Path(S3Path path){
        delegate = path
    }

    @Override
    void setTags(Map<String, String> tags) {
        this.tags = tags
    }

    @Override
    void setContentType(String type) {
        this.contentType = type
    }

    @Override
    void setStorageClass(String storageClass) {
        this.storageClass = storageClass
    }

    List<Tag> getTagsList() {
        // nothing found, just return
        if( tags==null )
            return Collections.emptyList()
        // create a list of Tag out of the Map
        List<Tag> result = new ArrayList<>()
        for( Map.Entry<String,String> entry : tags.entrySet()) {
            result.add( Tag.builder().key(entry.getKey()).value(entry.getValue()).build() )
        }
        return result;
    }

    String getContentType() {
        return contentType
    }

    String getStorageClass() {
        return storageClass
    }

    @Override
    FileSystem getFileSystem() {
        return delegate.getFileSystem()
    }

    @Override
    boolean isAbsolute() {
        return delegate.isAbsolute()
    }

    @Override
    Path getRoot() {
        return delegate.getRoot()
    }

    @Override
    Path getFileName() {
        return new NextflowS3Path(delegate.getFileName())
    }

    @Override
    Path getParent() {
        return new NextflowS3Path(delegate.getParent())
    }

    @Override
    int getNameCount() {
        return delegate.getNameCount()
    }

    @Override
    Path getName(int index) {
        return delegate.getName(index)
    }

    @Override
    Path subpath(int beginIndex, int endIndex) {
        return delegate.subpath(beginIndex,endIndex)
    }

    @Override
    boolean startsWith(Path other) {
        return delegate.startsWith( unwrapS3Path(other) )
    }

    @Override
    boolean startsWith(String other) {
        return delegate.startsWith(other)
    }

    @Override
    boolean endsWith(Path other) {
        return delegate.endsWith(unwrapS3Path(other))
    }

    @Override
    boolean endsWith(String other) {
        return delegate.endsWith(other)
    }

    @Override
    Path normalize() {
        return new NextflowS3Path(delegate.normalize())
    }

    @Override
    Path resolve(Path other) {
        return new NextflowS3Path(delegate.resolve(unwrapS3Path(other)))
    }

    @Override
    Path resolve(String other) {
        return new NextflowS3Path(delegate.resolve(other))
    }

    @Override
    Path resolveSibling(Path other) {
        return new NextflowS3Path(delegate.resolveSibling(unwrapS3Path(other)))
    }

    @Override
    Path resolveSibling(String other) {
        return new NextflowS3Path( delegate.resolveSibling(other))
    }

    @Override
    Path relativize(Path other) {
        return new NextflowS3Path( delegate.relativize( unwrapS3Path(other) ) )
    }

    @Override
    URI toUri() {
        return delegate.toUri()
    }

    @Override
    Path toAbsolutePath() {
        return new NextflowS3Path(delegate.toAbsolutePath())
    }

    @Override
    Path toRealPath(LinkOption... options) throws IOException {
        return delegate.toRealPath(options)
    }

    @Override
    File toFile() {
        return delegate.toFile()
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>[] events, WatchEvent.Modifier... modifiers) throws IOException {
        return delegate.register(watcher, events, modifiers)
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>... events) throws IOException {
        return delegate.register(watcher, events)
    }

    @Override
    Iterator<Path> iterator() {
        final iterator = delegate.iterator()
        return new Iterator<Path>() {
            @Override
            boolean hasNext() {
                return iterator.hasNext()
            }

            @Override
            Path next() {
                return new NextflowS3Path(iterator.next() as S3Path)
            }
        }
    }

    @Override
    int compareTo(Path other) {
        if( other instanceof NextflowS3Path )
            return delegate.compareTo((other as NextflowS3Path).toS3Path())
        return delegate.compareTo(other)
    }

    @Override
    String toString() {
        final bucket = delegate.bucketName()
        final key = delegate.getKey()
        return "/${bucket}/${key}"
    }

    @Override
    int hashCode() {
        return delegate.hashCode()
    }

    NextflowS3PathOpenOptions getOpenOptions(){
        return new NextflowS3PathOpenOptions(tagsList, storageClass, contentType)
    }

    S3Path toS3Path(){
        return delegate
    }

    static Path unwrapS3Path(Path other){
        if( other instanceof NextflowS3Path )
            return (other as NextflowS3Path).toS3Path()
        else
            return other
    }


}
