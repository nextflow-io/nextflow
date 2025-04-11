package nextflow.plugin

import com.github.tomakehurst.wiremock.junit.WireMockRule
import com.github.tomjankes.wiremock.WireMockGroovy
import dev.failsafe.FailsafeException
import org.junit.Rule
import org.pf4j.PluginRuntimeException
import spock.lang.Specification

import java.time.ZoneOffset
import java.time.ZonedDateTime

class HttpPluginRepositoryTest extends Specification {
    @Rule
    WireMockRule wiremock = new WireMockRule(0)

    def wm
    HttpPluginRepository unit

    def setup() {
        wm = new WireMockGroovy(wiremock.port())
        unit = new HttpPluginRepository("test-repo", new URI(wiremock.baseUrl()))
    }

    // ------------------------------------------------------------------------

    def 'prefetch metadata for plugin with no releases'() {
        given:
        wm.stub {
            request {
                method 'POST'
                url '/plugins/collect'
            }
            response {
                status 200
                body """{
                  "plugins": [
                    {
                      "id": "nf-fake",
                      "releases": []
                    }
                  ]
                }
                """
            }
        }

        when:
        unit.prefetch([new PluginSpec("nf-fake")])

        then:
        def plugins = unit.getPlugins()
        plugins.size() == 1
        def p1 = plugins.get("nf-fake")
        p1.id == "nf-fake"
        p1.releases.size() == 0
    }

    // ------------------------------------------------------------------------

    def 'prefetch plugin metadata with release'() {
        given:
        wm.stub {
            request {
                method 'POST'
                url '/plugins/collect'
            }
            response {
                status 200
                body """{
                  "plugins": [
                    {
                      "id": "nf-fake",
                      "releases": [
                        {
                          "version": "0.0.1",
                          "url": "http://example.com/fake-plugin/0.0.1/plugin.zip",
                          "date": "2025-03-21T10:40:36Z",
                          "sha512sum": "cbc4a7f0bc10c955ff10b85da85f4df8303e4174f1c44922c52a274a59afd786ee3bc685f332d1e0a7497d0b28b0cfc9939d052d6fa992607f43c870825f6caf",
                          "requires": ">=24.10.0"
                        }
                      ]
                    }
                  ]
                }
                """
            }
        }

        when:
        unit.prefetch([new PluginSpec("nf-fake")])

        then:
        def plugins = unit.getPlugins()
        plugins.size() == 1
        def p1 = plugins.get("nf-fake")
        p1.id == "nf-fake"
        p1.releases.size() == 1
        def r1 = p1.releases[0]
        r1.version == "0.0.1"
        r1.url == "http://example.com/fake-plugin/0.0.1/plugin.zip"
        r1.date == toDate(ZonedDateTime.of(2025, 3, 21, 10, 40, 36, 0, ZoneOffset.UTC))
        r1.sha512sum == "cbc4a7f0bc10c955ff10b85da85f4df8303e4174f1c44922c52a274a59afd786ee3bc685f332d1e0a7497d0b28b0cfc9939d052d6fa992607f43c870825f6caf"
        r1.requires == ">=24.10.0"
    }

    // ------------------------------------------------------------------------

    def 'handle prefetch error when metadata service is unavailable'() {
        setup:
        wiremock.stop()

        when:
        unit.prefetch([new PluginSpec("nf-fake")])

        then:
        def err = thrown FailsafeException
        err.message == "java.net.ConnectException: Failed to download plugins metadata"
    }

    // ------------------------------------------------------------------------

    def 'handle prefetch error when metadata service returns an error response'() {
        given:
        wm.stub {
            request {
                method 'POST'
                url '/plugins/collect'
            }
            response {
                status 500
                body "Server error!"
            }
        }

        when:
        unit.prefetch([new PluginSpec("nf-fake")])

        then:
        def err = thrown PluginRuntimeException
        err.message == "Failed to download plugin metadata: Server error!"
    }

    // ------------------------------------------------------------------------

    def 'handle prefetch error when metadata service sends back incorrectly formatted response'() {
        given:
        wm.stub {
            request {
                method: 'POST'
                url: '/plugins/collect'
            }
            response {
                status 200
                body """{
                  "not-plugins": [
                    {
                      "id": "nf-fake"
                    }
                  ]
                }
                """
            }
        }

        when:
        unit.prefetch([new PluginSpec("nf-fake")])

        then:
        def err = thrown PluginRuntimeException
        err.message == "Failed to download plugin metadata: Failed to parse response body"
    }

    // ------------------------------------------------------------------------

    def 'handle prefetch error caused by nextflow sending a bad request to metadata service'() {
        given:
        wm.stub {
            request {
                method: 'POST'
                url: '/plugins/collect'
            }
            response {
                status 400
                body """{
                  "type": "SOME_ERROR",
                  "message": "Unparseable request"
                }"""
            }
        }

        when:
        unit.prefetch([new PluginSpec("nf-fake")])

        then:
        def err = thrown PluginRuntimeException
        err.message == "Failed to download plugin metadata: SOME_ERROR - Unparseable request"
    }

    // ------------------------------------------------------------------------

    private static Date toDate(ZonedDateTime from) {
        return Date.from(from.toInstant())
    }
}
