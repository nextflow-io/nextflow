/*
 * Copyright 2018, University of Tübingen, Quantitative Biology Center (QBiC)
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
import groovy.util.logging.Slf4j

/**
 * Small and simple http client that sends POST requests
 * to a given URL. Currently used by the MessageObserver class
 * only.
 *
 * @author Sven Fillinger <sven.fillinger@qbic.uni-tuebingen.de>
 */
@Slf4j
class SimpleHttpClient {

    /**
     * Default user agent
     */
    private static String USER_AGENT = "Nextflow weblog"

    /**
     * Contains the HTTP POST target url
     */
    private String url

    /**
     * Contains the response code of a request
     */
    private int responseCode = -1

    /**
     * Contains the response string of a request
     */
    private String response

    /**
     * Setter for the target url
     * @param url
     */
    void setUrl(String url) {
        this.url = url
    }

    /**
     * Send a json formatted string as HTTP POST request
     * @param json Message content as JSON
     */
    void sendHttpMessage(String json) throws IllegalStateException, IllegalArgumentException{

        if (!this.url)
            throw new IllegalStateException("URL needs to be set!")

        if (!isJson(json)){
            throw new IllegalArgumentException("Message String needs to be a valid JSON object!")
        }

        // Open a connection to the target url
        def con = setUpConnection(this.url)
        // Make header settings
        con.setRequestMethod("POST")
        con.setRequestProperty("Content-Type", "application/json")
        con.setRequestProperty("User-Agent", USER_AGENT)
        con.setRequestProperty("Accept-Language", "en-US,en;q=0.5")
        con.setDoOutput(true)

        // Send POST request
        DataOutputStream output = new DataOutputStream(con.getOutputStream())
        output.writeBytes(json)
        output.flush()
        output.close()

        // Retrieve response code
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
     * Sets up the header for the request properly and creates
     * a connection to a target url
     * @param url The target url
     * @return A http connection
     */
    protected HttpURLConnection setUpConnection(String url){
        new URL(url).openConnection() as HttpURLConnection
    }

    /**
     * Getter for url
     * @return The current set url
     */
    String getUrl(){
        return this.url
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
        return this.userAgent
    }

    /**
     * Setter for the header user agent
     * @param agent A user agent for the HTTP request
     */
    void setUserAgent(String agent){
        this.userAgent = agent
    }



}
