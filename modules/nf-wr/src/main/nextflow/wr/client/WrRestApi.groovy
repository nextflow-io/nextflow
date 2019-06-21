/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.wr.client

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

import java.security.GeneralSecurityException
import java.security.KeyStore
import java.security.SecureRandom
import java.security.cert.Certificate
import java.security.cert.CertificateException
import java.security.cert.CertificateFactory
import java.security.cert.X509Certificate
import javax.net.ssl.HostnameVerifier
import javax.net.ssl.KeyManager
import javax.net.ssl.SSLContext
import javax.net.ssl.SSLSession
import javax.net.ssl.TrustManager
import javax.net.ssl.TrustManagerFactory
import javax.net.ssl.X509TrustManager
import javax.net.ssl.HttpsURLConnection
import groovy.json.JsonSlurper
import static groovy.json.JsonOutput.toJson

import nextflow.wr.executor.WrFileCopyStrategy

import nextflow.processor.TaskRun

/**
 * WrRestApi lets you interact with wr's REST API.
 *
 * See https://github.com/VertebrateResequencing/wr/wiki/REST-API
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 */
@Slf4j
@CompileStatic
class WrRestApi {
    
    static private String endpoint
    static private String token

    static private final String rest = "/rest/v1"
    static private final String jobs = "/jobs/"

    WrRestApi(String endpoint, String token, String cacertPath) {
        this.endpoint = endpoint
        this.token = token

        // have Groovy trust cacertPath
        def nullHostnameVerifier = [
            verify: { hostname, session -> true }
        ]
        HttpsURLConnection.setDefaultHostnameVerifier(nullHostnameVerifier as HostnameVerifier)

        def is = new File(cacertPath).newInputStream()

        TrustManager[] trustManagers = null;
        char[] password = null; // any password will work.
        KeyStore caKeyStore = KeyStore.getInstance(KeyStore.getDefaultType());
        caKeyStore.load(null, password);
        CertificateFactory certificateFactory = CertificateFactory.getInstance("X.509");
        Collection<? extends Certificate> certificates = certificateFactory.generateCertificates(is);
        if (certificates.isEmpty()) {
            throw new IllegalArgumentException("expected non-empty set of trusted certificates");
        }
        int index = 0;
        certificates.each {
            String certificateAlias = "ca" + Integer.toString(index++);
            caKeyStore.setCertificateEntry(certificateAlias, it);
        }
        TrustManagerFactory trustManagerFactory = TrustManagerFactory.getInstance(TrustManagerFactory.getDefaultAlgorithm());
        trustManagerFactory.init(caKeyStore);
        trustManagers = trustManagerFactory.getTrustManagers();

        SSLContext sslContext = SSLContext.getInstance("SSL");
        sslContext.init(null, trustManagers, new SecureRandom());
        HttpsURLConnection.setDefaultSSLSocketFactory(sslContext.getSocketFactory());

        is.close()
    }

    /**
     * Adds new commands to the job queue. Takes a List of Lists, where the
     * sub lists are [String cmd, TaskRun task, WrFileCopyStrategy copyStrategy]
     * Returns the jobs.
     */
    List<Map> add(List<List> argList) {
        // curl --cacert ~/.wr_production/ca.pem -H "Content-Type: application/json" -H "Authorization: Bearer [token]" -X POST -d '[{"cmd":"sleep 5 && echo mymsg && false","memory":"5M","cpus":1},{"cmd":"sleep 5","cpus":1}]' 'https://localhost:11302/rest/v1/jobs/?rep_grp=myid&cpus=2&memory=3G&time=5s'

        List<Map> jsonArgs = []

        argList.each {
            String cmd = it[0]
            TaskRun task = it[1]
            WrFileCopyStrategy copyStrategy = it[2]

            String cwd = ""
            boolean cwdMatters = true
            switch (copyStrategy.getWorkDirScheme()) {
                case "file":
                    cwd = task.getWorkDirStr()
                    break
                case "s3":
                    cwd = "/tmp"
                    cwdMatters = false
                    break
                default:
                    throw new IllegalArgumentException("Unsupported scheme for work dir: ${task.getWorkDirStr()}")
            }

            String grp = "[nextflow] $task.processor.name"

            // calculate mounts option
            def reads = copyStrategy.inputBuckets()
            List<Map> m = []
            for (path in reads) {
                m << ["Mount":copyStrategy.getInputMountPath() + "/" + path, "Targets":[["Path":path]]]
            }
            def write = copyStrategy.outputBucket()
            if (write) {
                m << ["Mount":copyStrategy.getOutputMountPath(),"Targets":[["Path":write,"Write":true]]]
            }

            Integer cpus = 1
            if ( task.config.cpus && task.config.cpus > 0 ) {
                cpus = task.config.cpus
            }

            Integer override = 0

            String mem = ""
            if( task.config.getMemory() ) {
                mem = String.valueOf(task.config.getMemory().toUnit('MB')) + "M"
                override = 2
            }

            String t = ""
            if( task.config.time ) {
                t = task.config.getTime().format('Hh')
                override = 2
            }

            Long d = 0
            if( task.config.getDisk() ) {
                def disk = task.config.getDisk()
                d = disk.toUnit('GB')
                override = 2
            }

            List<String> limits = []
            if ( task.config.maxForks ) {
                limits << sprintf('%s:%d', task.processor.name, task.config.maxForks)
            }

            List<Map> behaviours = []
            behaviours << ["cleanup":true]

            // *** what about cloud opts like image and flavor?
            // Turn on docker monitoring if docker container is being used?
            // Does BashWrapperBuilder handle everything to do with env vars?

            Map args = [cmd: cmd, cwd: cwd, cwd_matters: cwdMatters, rep_grp: grp, req_grp: grp, limit_grps: limits, override: override, retries: 0, rerun: true, cpus: cpus, memory: mem, time: t, disk: d, mounts: m, on_exit: behaviours]
            // log.debug "[wr] add args: $args"
            jsonArgs << args
        }

        def response = postJson(jobs, jsonArgs)
        // log.debug("made a POST request")
        return parseJobsFromJson(response)
    }

    /**
     * Get the status of a previously add()ed job(s). Returns a List of Map with
     * details of each job, including Key, State, Exited and Exitcode.
     */
    List<Map> status(String ids) {
        // curl --cacert ~/.wr_production/ca.pem -H "Authorization: Bearer [token]" 'https://localhost:11302/rest/v1/jobs/58cef10e7a340c3b7fa09ea304a3cb98'

        // log.debug "working on ids $ids"

        def get = authenticatedConnection("$jobs/$ids")
        get.setRequestMethod("GET")
        def getRC = get.getResponseCode()
        // log.debug("made a GET request")
        if (!getRC.equals(200)) {
            log.debug("[wr] failed to get command details from wr manager")
            return null
        }
        return parseJobsFromJson(get.getInputStream().getText())
    }

    /**
     * Cancel the given job. If it's running it will be killed. If it's lost
     * it will be buried. If it's otherwise incomplete, it will be removed from
     * wr's queue.
     */
    void cancel(String id) {
        // like status, but send DELETE instead of GET. Also state parameter
        // required to be running|lost|deletable

        // *** is there a way to avoid needing to get the current state of the
        // job first? Does our caller only call us for running jobs?
        List jobs = status(id)
        Map job = jobs[0]
        String state
        if (job.State == 'running') {
            state = 'running'
        } else if (job.State == 'lost') {
            state = 'lost'
        } else if (job.State != 'complete') {
            state = 'deletable'
        }

        def delete = authenticatedConnection("$jobs/$id?state=$state")
        delete.setRequestMethod("DELETE")
        def deleteRC = delete.getResponseCode()
        if (!deleteRC.equals(200) || !deleteRC.equals(202)) {
            log.debug("[wr] failed to cancel command held in wr manager")
        }
    }

    private HttpsURLConnection authenticatedConnection(String leaf) {
        String url = "$endpoint$rest$leaf"
        // log.debug("[wr] url: $url")
        def post = new URL(url).openConnection() as HttpsURLConnection
        post.setRequestProperty("Content-Type", "application/json")
        post.setRequestProperty("Authorization", "Bearer $token")
        return post
    }

    private HttpsURLConnection postConnection(String leaf) {
        def post = authenticatedConnection(leaf)
        post.setRequestMethod("POST")
        post.setDoOutput(true)
        return post
    }

    private String postJson(String leaf, List<Map> args) {
        def json = toJson(args)
        // log.debug("[wr] posting JSON: $json")

        def post = postConnection(leaf)
        post.getOutputStream().write(json.getBytes("UTF-8"))
        def postRC = post.getResponseCode()
        // log.debug(postRC as String)
        if (!postRC.equals(201)) {
            log.debug("[wr] failed to send new command details to wr manager")
            return null
        }
        return post.getInputStream().getText()
    }

    private List<Map> parseJobsFromJson(String json) {
        List<Map> jobs = new JsonSlurper().parseText(json) as List<Map>
        // log.debug "got ${jobs.size()} jobs"
        return jobs
    }

    private String parseIdFromJson(String json) {
        return extractFirstJob(parseJobsFromJson(json)).Key
    }

    private Map extractFirstJob(List<Map> jobs) {
        if (jobs.size() != 1) {
            log.debug("[wr] expected 1 job, got ${jobs.size()}")
        }
        Map job = jobs[0]
        return job
    }

}