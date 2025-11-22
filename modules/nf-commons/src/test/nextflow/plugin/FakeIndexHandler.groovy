package nextflow.plugin

import com.sun.net.httpserver.Headers
import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpHandler

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FakeIndexHandler implements HttpHandler {

    String PLUGINS = '''
        [
          {
            "id": "nf-console",
            "releases": [
              {
                "version": "1.0.0",
                "url": "http://localhost:9900/download/1.0.0/nf-console-1.0.0.zip",
                "date": "2021-01-04T17:49:04.62+01:00",
                "sha512sum": "5875f1fd4c24ab9ecded06fa43a5b87b4e45f70085fe5c8f740da66b318ef3b45005b3dd9a15e69fff20bde52b4d42db8d25163ac49878f64d784d717102f295"
              }
            ]
          },
          {
            "id": "nf-hello",
            "releases": [
              {
                "version": "0.2.0",
                "url": "http://localhost:9900/download/1.0.0/nf-hello-0.2.0.zip",
                "date": "2021-01-04T17:49:04.62+01:00",
                "sha512sum": "c29fb06f785becc15117fc0587d8f50a313b92f07604a4770f1db9bf55f67a2c544a70588cf93fc1738cdaa43bb083ef9dfbf60e8dc1a7c2fdbce47cbb474441"
              }
            ]
          },
          {
            "id": "nf-amazon",
            "releases": [
              {
                "version": "1.0.0",
                "url": "https://github.com/nextflow-io/nf-amazon/releases/download/1.0.0/nf-amazon-1.0.0.zip",
                "date": "2021-01-04T17:49:04.62+01:00",
                "sha512sum": "5abe4cbc643ca0333cba545846494b17488d19d17049bf9d4d9ac9ec108352b54f147abb15c48581b44fdd09097e3607aed77e2eae91c07da32c2cd71f4f852f"
              }
            ]
          },
          {
            "id": "nf-google",
            "releases": [
              {
                "version": "1.0.0",
                "url": "https://github.com/nextflow-io/nf-google/releases/download/1.0.0/nf-google-1.0.0.zip",
                "date": "2021-01-04T17:49:04.62+01:00",
                "sha512sum": "023037d6daedc5e61b28bbd881843b41dffe5c5a81c982569288055e60568bba55a8704cdc237669dcc7b148af5c30674e74338cb48da7ab5674db43e3591137"
              }
            ]
          },
          {
            "id": "nf-tower",
            "releases": [
              {
                "version": "1.0.0",
                "url": "https://github.com/nextflow-io/nf-tower/releases/download/1.0.0/nf-tower-1.0.0.zip",
                "date": "2021-01-04T17:49:04.62+01:00",
                "sha512sum": "c00435164b4f16e548df4bea6f0e8c271dd8d0eaf28b10f29758d0a93aca5c426507c913a2ed0828e1a79e2c01d8f06644df871ccf0de2224f56745110f4a923"
              }
            ]
          }
        ]
'''

    File zip = new File('src/testResources/nf-console-1.0.0.zip')
    File zip2 = new File('src/testResources/nf-hello-0.2.0.zip')

    @Override
    void handle(HttpExchange request) throws IOException {

        def path = request.requestURI.toString()
        if( path.endsWith('plugins.json') ) {
            replyWithJson(request, PLUGINS)
        }
        else if( path.endsWith(zip.name)) {
            replyWithZip(request, zip)
        }
        else if( path.endsWith(zip2.name)) {
            replyWithZip(request, zip2)
        }
        else {
            new IllegalArgumentException("Invalid request: $path")
        }
    }

    void replyWithJson(HttpExchange request, String json) {
        Headers header = request.getResponseHeaders()
        header.set("Content-Type", "text/plain")
        request.sendResponseHeaders(200, json.size())

        OutputStream os = request.getResponseBody();
        os.write(json.bytes);
        os.close();
    }

    void replyWithZip(HttpExchange request, File file) {
        Headers header = request.getResponseHeaders()
        header.set("Content-Type", "text/plain")
        request.sendResponseHeaders(200, file.size())

        OutputStream os = request.getResponseBody();
        os.write(file.bytes);
        os.close();
    }
}
