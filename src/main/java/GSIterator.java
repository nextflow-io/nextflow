import com.google.common.base.Preconditions;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class GSIterator implements Iterator<Path> {

    private GSFileSystem gsFileSystem;
    private String bucket;
    private String key;

    private Iterator<GSPath> it;

    public GSIterator(GSFileSystem gsFileSystem, String bucket, String key) {

        Preconditions.checkArgument(key != null && key.endsWith("/"), "key %s should be ended with slash '/'", key);

        this.gsFileSystem = gsFileSystem;
        this.bucket = bucket;
        this.key = key.length() == 1 ? "" : key;
    }

    @Override
    public void remove() { throw new UnsupportedOperationException(); }

    @Override
    public GSPath next() { return getIterator().next(); }

    @Override
    public boolean hasNext() { return getIterator().hasNext(); }

    private Iterator<GSPath> getIterator() {
        if (it == null) {
            List<GSPath> listPath = new ArrayList<>();

            // iterate over this list
            ObjectListing
        }
    }
}
