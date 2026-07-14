/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package test

import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * A minimal raw-socket forward HTTP proxy requiring Basic authentication,
 * for testing authenticated proxy support.
 *
 * Requests without a valid {@code Proxy-Authorization} header are rejected
 * with a {@code 407} response. Authenticated plain HTTP requests (absolute-form
 * request target) are answered directly by the proxy with the configured canned
 * response, without contacting any origin server. Authenticated {@code CONNECT}
 * requests are tunnelled to the requested target, allowing HTTPS traffic to
 * flow through the proxy.
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@Slf4j
@CompileStatic
class MockAuthProxyServer implements Closeable {

    private final String username
    private final String password
    private final String expectedAuth
    private ServerSocket server
    private Thread acceptor
    private volatile boolean running

    /** canned response returned to authenticated plain HTTP requests */
    String responseBody = 'OK'
    String responseContentType = 'text/plain'

    final AtomicInteger requestCount = new AtomicInteger()
    final AtomicInteger unauthorizedCount = new AtomicInteger()
    final AtomicInteger proxiedCount = new AtomicInteger()
    final AtomicInteger connectCount = new AtomicInteger()

    MockAuthProxyServer(String username, String password) {
        this.username = username
        this.password = password
        this.expectedAuth = 'Basic ' + Base64.getEncoder().encodeToString("${username}:${password}".getBytes('UTF-8'))
    }

    /**
     * Environment variables routing all HTTP traffic through this proxy with the
     * given credentials - by default the ones expected by the proxy. Meant to be
     * applied with {@code SysEnv.push()}, which replaces the environment map as
     * seen by the code under test, so no other proxy settings can be inherited
     * from the host environment
     */
    Map<String,String> proxyEnv(String user=username, String password=this.password) {
        final result = new HashMap<String,String>()
        result.put('HTTP_PROXY', "http://${user}:${password}@${host}:${port}".toString())
        result.put('HTTPS_PROXY', result.get('HTTP_PROXY'))
        return result
    }

    MockAuthProxyServer start() {
        server = new ServerSocket(0, 50, InetAddress.getByName('127.0.0.1'))
        running = true
        acceptor = Thread.startDaemon("MockAuthProxyServer-acceptor") {
            while( running ) {
                try {
                    final socket = server.accept()
                    Thread.startDaemon("MockAuthProxyServer-worker") { handle(socket) }
                }
                catch( IOException e ) {
                    if( running )
                        log.warn "Mock proxy accept error", e
                }
            }
        }
        return this
    }

    int getPort() { server.getLocalPort() }

    String getHost() { '127.0.0.1' }

    @Override
    void close() {
        running = false
        server?.close()
    }

    private void handle(Socket socket) {
        try {
            socket.withCloseable {
                final input = socket.getInputStream()
                final output = socket.getOutputStream()
                final requestLine = readLine(input)
                if( !requestLine )
                    return
                final headers = new LinkedHashMap<String,String>()
                String line
                while( (line=readLine(input)) ) {
                    final p = line.indexOf(':')
                    if( p!=-1 )
                        headers.put(line.substring(0,p).trim().toLowerCase(), line.substring(p+1).trim())
                }
                requestCount.incrementAndGet()
                log.debug "Mock proxy request: $requestLine"

                if( !isAuthorized(headers) ) {
                    unauthorizedCount.incrementAndGet()
                    final resp = "HTTP/1.1 407 Proxy Authentication Required\r\n" +
                            "Proxy-Authenticate: Basic realm=\"test\"\r\n" +
                            "Content-Length: 0\r\n" +
                            "Connection: close\r\n\r\n"
                    output.write(resp.getBytes('ISO-8859-1'))
                    output.flush()
                    return
                }

                if( requestLine.startsWith('CONNECT ') ) {
                    connectCount.incrementAndGet()
                    tunnel(requestLine, socket, output)
                }
                else {
                    proxiedCount.incrementAndGet()
                    final body = responseBody.getBytes('UTF-8')
                    final resp = "HTTP/1.1 200 OK\r\n" +
                            "Content-Type: ${responseContentType}\r\n" +
                            "Content-Length: ${body.length}\r\n" +
                            "Connection: close\r\n\r\n"
                    output.write(resp.getBytes('ISO-8859-1'))
                    output.write(body)
                    output.flush()
                }
            }
        }
        catch( IOException e ) {
            log.debug "Mock proxy connection error - ${e.message}"
        }
    }

    private boolean isAuthorized(Map<String,String> headers) {
        return headers.get('proxy-authorization') == expectedAuth
    }

    private void tunnel(String requestLine, Socket socket, OutputStream output) {
        // request line is expected as `CONNECT host:port HTTP/1.1` - the connection
        // is always tunnelled to the loopback address since the target host names
        // used by the tests are synthetic and cannot be resolved
        final target = requestLine.tokenize(' ')[1]
        final port = target.substring(target.lastIndexOf(':')+1) as int
        new Socket('127.0.0.1', port).withCloseable { upstream ->
            output.write("HTTP/1.1 200 Connection established\r\n\r\n".getBytes('ISO-8859-1'))
            output.flush()
            final t1 = Thread.startDaemon("MockAuthProxyServer-up") { pump(socket.getInputStream(), upstream.getOutputStream()) }
            pump(upstream.getInputStream(), socket.getOutputStream())
            t1.join(5_000)
        }
    }

    private static void pump(InputStream input, OutputStream output) {
        try {
            final buffer = new byte[8192]
            int n
            while( (n=input.read(buffer)) != -1 ) {
                output.write(buffer, 0, n)
                output.flush()
            }
        }
        catch( IOException e ) {
            // ignore - connection closed by either side
        }
    }

    private static String readLine(InputStream input) {
        final buffer = new ByteArrayOutputStream()
        int ch
        while( (ch=input.read()) != -1 ) {
            if( ch == '\n' )
                break
            if( ch != '\r' )
                buffer.write(ch)
        }
        return buffer.toString('ISO-8859-1')
    }
}
