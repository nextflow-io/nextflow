/*
 * Copyright 2018, University of Tübingen, Quantitative Biology Center (QBiC)
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.util

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const

/**
 * Small and simple http client that sends POST requests
 * to a given URL. Currently used by the MessageObserver class
 * only.
 *
 * @author
 *  Sven Fillinger <sven.fillinger@qbic.uni-tuebingen.de>
 *  Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SimpleHttpClient {

    public static int DEFAULT_BACK_OFF_BASE = 3
    public static int DEFAULT_BACK_OFF_DELAY = 250

    /**
     * Default user agent
     */
    private static String DEF_USER_AGENT = "Nextflow/$Const.APP_VER"

    /**
     * Contains the response code of a request
     */
    private int responseCode = -1

    /**
     * Contains the response string of a request
     */
    private String response

    /**
     * Basic auth token
     */
    private String authToken

    /**
     * Http user agent
     */
    private String userAgent = DEF_USER_AGENT

    private int errorCount

    private int maxRetries

    private int backOffBase = DEFAULT_BACK_OFF_BASE

    private int backOffDelay = DEFAULT_BACK_OFF_DELAY

    /**
     * Send a json formatted string as HTTP POST request
     * @param json Message content as JSON
     */
    void sendHttpMessage(String url, String json, String method = 'POST') throws IllegalStateException, IllegalArgumentException{

        if (!url)
            throw new IllegalStateException("URL needs to be set!")

        // reset the error count 
        errorCount = 0

        while( true ) {
            // Open a connection to the target url
            def con = getHttpConnection(url)
            // Make header settings
            con.setRequestMethod(method)
            con.setRequestProperty("Content-Type", "application/json")
            con.setRequestProperty("User-Agent", userAgent)
            if( authToken )
                con.setRequestProperty("Authorization","Basic ${authToken.bytes.encodeBase64()}")

            con.setDoOutput(true)

            // Send POST request
            if( json ) {
                DataOutputStream output = new DataOutputStream(con.getOutputStream())
                output.writeBytes(json)
                output.flush()
                output.close()
            }

            // Retrieve response code
            try {
                this.responseCode = con.getResponseCode()

                // Retrieve response message
                BufferedReader reader = new BufferedReader(new InputStreamReader(con.getInputStream()))
                String lineContent
                StringBuffer tmpResponse = new StringBuffer()

                while ((lineContent = reader.readLine()) != null){
                    tmpResponse.append(lineContent)
                }

                reader.close()

                this.response = tmpResponse.toString()
                break
            }
            catch( IOException e ) {
                // capture error code
                // https://stackoverflow.com/a/18462721/395921
                this.responseCode = con.getResponseCode()
                this.response = (con.getErrorStream() ?: con.getInputStream())?.text
                this.errorCount +=1
                if( responseCode < 500 || this.errorCount > maxRetries )
                    throw e

                final delay = (Math.pow(backOffBase, errorCount) as long) * backOffDelay
                log.debug "Got HTTP error=$responseCode waiting for ${delay}ms (errorCount=$errorCount)"
                Thread.sleep(delay)
            }
        }

    }

    protected HttpURLConnection getHttpConnection(String url) {
        new URL(url).openConnection() as HttpURLConnection
    }

    /**
     * Validate JSON string
     * @param s Putative JSON object
     * @return true if JSON, false else
     */
    protected boolean isJson(String s){
        try{
            new JsonSlurper().parseText(s)
        } catch (Exception e){
            log.error(e.message, e)
            return false
        }
        return true
    }

    /**
     * Getter for the response code
     * @return The HTTP response code
     */
    int getResponseCode(){
        return this.responseCode
    }

    /**
     * Getter for the response message
     * @return The HTTP response message
     */
    String getResponse(){
        return this.response
    }

    /**
     * Getter for the current user agent
     * @return The HTTP header user agent
     */
    String getUserAgent(){
        return userAgent
    }

    /**
     * Setter for the header user agent
     * @param agent A user agent for the HTTP request
     */
    SimpleHttpClient setUserAgent(String agent){
        userAgent = agent
        return this
    }

    SimpleHttpClient setAuthToken(String tkn) {
        authToken = tkn
        return this
    }

    SimpleHttpClient setMaxRetries(int value) {
        this.maxRetries = value
        return this
    }

    SimpleHttpClient setBackOffBase(int value ) {
        this.backOffBase = value
        return this
    }

    SimpleHttpClient setBackOffDelay(int value ) {
        this.backOffDelay= value
        return this
    }

}
