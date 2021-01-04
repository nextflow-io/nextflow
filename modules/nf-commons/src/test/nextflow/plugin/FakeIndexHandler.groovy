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
                "sha512sum": "4e681a539a0a902c1d93bce04cb2a5e49ed286581bc1ea0202044d040566f583aa0efabe0865d351eb284e4962020bc23d62075f52b04a06ab655b53e8aeea34"
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
                "sha512sum": "4e681a539a0a902c1d93bce04cb2a5e49ed286581bc1ea0202044d040566f583aa0efabe0865d351eb284e4962020bc23d62075f52b04a06ab655b53e8aeea34"
              }
            ]
          }
        ]
'''

    File zip = new File('src/testResources/nf-console-1.0.0.zip')

    @Override
    void handle(HttpExchange request) throws IOException {

        def path = request.requestURI.toString()
        if( path.endsWith('plugins.json') ) {
            replyWithJson(request, PLUGINS)
        }
        else if( path.endsWith(zip.name)) {
            replyWithZip(request, zip)
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
