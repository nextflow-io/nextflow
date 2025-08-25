package nextflow.plugin


import java.time.ZoneOffset
import java.time.ZonedDateTime

import com.github.tomakehurst.wiremock.junit.WireMockRule
import com.github.tomjankes.wiremock.WireMockGroovy
import nextflow.BuildInfo
import org.junit.Rule
import org.pf4j.PluginRuntimeException
import spock.lang.Specification

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
                method 'GET'
                url "/v1/plugins/dependencies?plugins=nf-fake&nextflowVersion=${BuildInfo.version}"
            }
            response {
                status 200
                body """{
                  "plugins": [
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
        def plugins = unit.getPlugins()
        plugins.size() == 0
    }

    // ------------------------------------------------------------------------

    def 'prefetch plugin metadata with release'() {
        given:
        wm.stub {
            request {
                method 'GET'
                url "/v1/plugins/dependencies?plugins=nf-fake&nextflowVersion=${BuildInfo.version}"
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
        def err = thrown PluginRuntimeException
        err.message.startsWith("Unable to connect to http://localhost")
    }

    // ------------------------------------------------------------------------

    def 'handle prefetch error when metadata service returns an error response'() {
        given:
        wm.stub {
            request {
                method 'GET'
                url "/v1/plugins/dependencies?plugins=nf-fake&nextflowVersion=${BuildInfo.version}"
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
        err.message == """\
            Invalid response while fetching plugin metadata from: http://localhost:${wiremock.port()}/v1/plugins/dependencies?plugins=nf-fake&nextflowVersion=${BuildInfo.version}
            - http status: 500
            - response   : Server error!""".stripIndent()
    }

    // ------------------------------------------------------------------------

    def 'handle prefetch error when metadata service sends back incorrectly formatted response'() {
        given:
        wm.stub {
            request {
                method 'GET'
                url "/v1/plugins/dependencies?plugins=nf-fake&nextflowVersion=${BuildInfo.version}"
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
        err.message.startsWith("Unexpected error while fetching plugin metadata from: http://localhost")
    }

    // ------------------------------------------------------------------------

    def 'handle prefetch error caused by nextflow sending a bad request to metadata service'() {
        given:
        wm.stub {
            request {
                method 'GET'
                url "/v1/plugins/dependencies?plugins=nf-fake&nextflowVersion=${BuildInfo.version}"
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
        err.message.startsWith("Invalid response while fetching plugin metadata from: http://localhost")
    }

    // ------------------------------------------------------------------------

    def 'test mapToPluginInfo with comprehensive date conversion validation'() {
        given:
        // Test the mapToPluginInfo method with real date conversion using GsonEncoder
        def utcDateStr = "2023-12-25T14:30:45Z"
        def estDateStr = "2023-06-15T09:15:30-05:00"
        def cestDateStr = "2023-08-10T16:45:00+02:00"
        
        wm.stub {
            request {
                method 'GET'
                url "/v1/plugins/dependencies?plugins=date-test-plugin&nextflowVersion=${BuildInfo.version}"
            }
            response {
                status 200
                body """{
                  "plugins": [
                    {
                      "id": "date-test-plugin",
                      "projectUrl": "https://example.com/test-plugin",
                      "provider": "Test Provider",
                      "releases": [
                        {
                          "version": "1.0.0",
                          "url": "https://example.com/plugin-1.0.0.zip",
                          "date": "${utcDateStr}",
                          "sha512sum": "hash1",
                          "requires": ">=20.0.0"
                        },
                        {
                          "version": "1.1.0", 
                          "url": "https://example.com/plugin-1.1.0.zip",
                          "date": "${estDateStr}",
                          "sha512sum": "hash2",
                          "requires": ">=21.0.0"
                        },
                        {
                          "version": "1.2.0",
                          "url": "https://example.com/plugin-1.2.0.zip", 
                          "date": "${cestDateStr}",
                          "sha512sum": "hash3",
                          "requires": ">=22.0.0"
                        },
                        {
                          "version": "1.3.0",
                          "url": "https://example.com/plugin-1.3.0.zip",
                          "sha512sum": "hash4",
                          "requires": ">=23.0.0"
                        }
                      ]
                    }
                  ]
                }
                """
            }
        }

        when:
        unit.prefetch([new PluginSpec("date-test-plugin")])

        then:
        def plugins = unit.getPlugins()
        plugins.size() == 1
        
        def pluginInfo = plugins.get("date-test-plugin")
        
        // Verify basic plugin mapping through the mapToPluginInfo method
        pluginInfo.id == "date-test-plugin"
        pluginInfo.projectUrl == "https://example.com/test-plugin"
        pluginInfo.provider == "Test Provider"
        pluginInfo.releases.size() == 4
        pluginInfo.releases != null // Verify never null
        
        // Verify UTC date conversion (Z suffix)
        def release1 = pluginInfo.releases[0]
        release1.version == "1.0.0"
        release1.date == toDate(ZonedDateTime.of(2023, 12, 25, 14, 30, 45, 0, ZoneOffset.UTC))
        release1.sha512sum == "hash1"
        release1.requires == ">=20.0.0"
        
        // Verify EST date conversion (-05:00 offset)
        def release2 = pluginInfo.releases[1]
        release2.version == "1.1.0"
        // 09:15 EST (-5 hours) = 14:15 UTC
        release2.date == toDate(ZonedDateTime.of(2023, 6, 15, 14, 15, 30, 0, ZoneOffset.UTC))
        release2.sha512sum == "hash2"
        release2.requires == ">=21.0.0"
        
        // Verify CEST date conversion (+02:00 offset)
        def release3 = pluginInfo.releases[2]
        release3.version == "1.2.0"
        // 16:45 CEST (+2 hours) = 14:45 UTC
        release3.date == toDate(ZonedDateTime.of(2023, 8, 10, 14, 45, 0, 0, ZoneOffset.UTC))
        release3.sha512sum == "hash3"
        release3.requires == ">=22.0.0"
        
        // Verify null date handling (missing date field)
        def release4 = pluginInfo.releases[3]
        release4.version == "1.3.0"
        release4.date == null
        release4.sha512sum == "hash4"
        release4.requires == ">=23.0.0"
    }


    // ------------------------------------------------------------------------

    private static Date toDate(ZonedDateTime from) {
        return Date.from(from.toInstant())
    }
}
