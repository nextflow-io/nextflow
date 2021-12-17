package com.upplication.s3fs.ng;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.model.GetObjectRequest;
import com.amazonaws.services.s3.model.S3Object;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;

public class DownloadOnDemandInputStream extends InputStream {

    private static final Logger log = LoggerFactory.getLogger(DownloadOnDemandInputStream.class);

    private final AmazonS3 s3;
    private final GetObjectRequest req;

    private InputStream parent;

    public DownloadOnDemandInputStream(GetObjectRequest req, AmazonS3 s3) {
        this.req = req;
        this.s3 = s3;
    }

    private InputStream parentInputStream() throws IOException {
        if (this.parent == null) {
            try (S3Object chunk = s3.getObject(req)) {
                final long[] range = req.getRange();
                final String chunkRef = range == null ?
                        String.format("part=%s", req.getPartNumber()) :
                        String.format("range=%s..%s", range[0], range[1]);
                final String path = "s3://" + req.getBucketName() + '/' + req.getKey();
                log.trace("Download chunk {}; path={}", chunkRef, path);
                try {
                    this.parent = chunk.getObjectContent();
                } catch (Throwable e) {
                    String msg = String.format("Failed to download chunk %s; path=%s", chunkRef, path);
                    throw new IOException(msg, e);
                }
            }
        }
        return this.parent;
    }

    @Override
    public int read() throws IOException {
        return parentInputStream().read();
    }

    @Override
    public int read(byte b[], int off, int len) throws IOException {
        return parentInputStream().read(b, off, len);
    }
}
