// Copyright (C) 2013 DNAnexus, Inc.
//
// This file is part of dx-toolkit (DNAnexus platform client libraries).
//
//   Licensed under the Apache License, Version 2.0 (the "License"); you may
//   not use this file except in compliance with the License. You may obtain a
//   copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
//   WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
//   License for the specific language governing permissions and limitations
//   under the License.

package com.dnanexus;

import com.fasterxml.jackson.core.*;
import com.fasterxml.jackson.databind.*;
import org.apache.http.*;
import org.apache.http.client.*;
import org.apache.http.util.*;
import org.apache.http.entity.*;
import org.apache.http.client.methods.*;
import org.apache.http.impl.client.DefaultHttpClient;
import org.apache.http.client.ClientProtocolException;
import java.io.*;
import org.apache.commons.io.IOUtils;

/**
 * Class for making a raw DNAnexus API call via HTTP.
 */
public class DXHTTPRequest {
    private final JsonNode securityContext;
    private final String apiserver;
    private final DefaultHttpClient httpclient;

    private static int NUM_RETRIES = 5;

    private static DXEnvironment env = DXEnvironment.create();

    public DXHTTPRequest() throws Exception {
        this.securityContext = env.getSecurityContext();
        this.apiserver = env.getApiserverPath();
        this.httpclient = new DefaultHttpClient();
    }

    private String errorMessage(String method, String resource, String errorString,
                                int retryWait, int nextRetryNum, int maxRetries) {
        String baseError = method + " " + resource + ": " + errorString + ".";
        if (nextRetryNum <= maxRetries) {
            return baseError + "  Waiting " + retryWait + " seconds before retry " + nextRetryNum
                + " of " + maxRetries;
        } else {
            return baseError;
        }
    }

    /**
     * Holds either the raw text of a response or a parsed JSON version of it.
     */
    private static class ParsedResponse {
        public final String responseText;
        public final JsonNode responseJson;

        public ParsedResponse(String responseText, JsonNode responseJson) {
            this.responseText = responseText;
            this.responseJson = responseJson;
        }
    }

    /**
     * Issues a request against the specified resource and returns either the
     * text of the response or the parsed JSON of the response (depending on
     * whether parseResponse is set).
     */
    private ParsedResponse requestImpl(String resource, String data, boolean parseResponse) throws Exception {
        HttpPost request = new HttpPost(apiserver + resource);

        request.setHeader("Content-Type", "application/json");
        request.setHeader("Authorization", securityContext.get("auth_token_type").textValue()
                          + " " + securityContext.get("auth_token").textValue());
        request.setEntity(new StringEntity(data));

        // Retry with exponential backoff
        int timeout = 1;

        for (int i = 0; i <= NUM_RETRIES; i++) {
            HttpResponse response = null;
            boolean okToRetry = false;

            try {
                response = httpclient.execute(request);
            } catch (ClientProtocolException e) {
                System.err.println(errorMessage("POST", resource, e.toString(), timeout, i + 1,
                                                NUM_RETRIES));
            } catch (IOException e) {
                System.err.println(errorMessage("POST", resource, e.toString(), timeout, i + 1,
                                                NUM_RETRIES));
            }

            if (response != null) {
                int statusCode = response.getStatusLine().getStatusCode();

                HttpEntity entity = response.getEntity();

                if (statusCode == HttpStatus.SC_OK) {
                    // 200 OK

                    byte[] value = EntityUtils.toByteArray(entity);
                    int realLength = value.length;
                    if (entity.getContentLength() >= 0 && realLength != entity.getContentLength()) {
                        String errorStr = "Received response of " + realLength
                            + " bytes but Content-Length was " + entity.getContentLength();
                        System.err.println(errorMessage("POST", resource, errorStr, timeout, i + 1,
                                                        NUM_RETRIES));
                    } else {
                        if (parseResponse) {
                            JsonNode responseJson = null;
                            try {
                                responseJson = DXJSON.parseJson(new String(value, "UTF-8"));
                            } catch (JsonProcessingException e) {
                                if (entity.getContentLength() < 0) {
                                    // content-length was not provided, and the
                                    // JSON could not be parsed. Retry since
                                    // this is a streaming request from the
                                    // server that probably just encountered a
                                    // transient error.
                                } else {
                                    throw e;
                                }
                            }
                            if (responseJson != null) {
                                return new ParsedResponse(null, responseJson);
                            }
                        } else {
                            return new ParsedResponse(new String(value, "UTF-8"), null);
                        }
                    }
                } else {
                    // Non-200 status codes.

                    // 500 InternalError should get retried. 4xx errors should
                    // be considered not recoverable.
                    if (statusCode < 500) {
                        throw new Exception(EntityUtils.toString(entity));
                    } else {
                        System.err.println(errorMessage("POST", resource, EntityUtils.toString(entity),
                                                        timeout, i + 1, NUM_RETRIES));
                    }
                }
            }

            if (i < NUM_RETRIES) {
                Thread.sleep(timeout * 1000);
                timeout *= 2;
            }
        }

        throw new Exception("POST " + resource + " failed");
    }

    /**
     * Issues a request against the specified resource and returns the result
     * as a String.
     */
    public String request(String resource, String data) throws Exception {
        return requestImpl(resource, data, false).responseText;
    }

    /**
     * Issues a request against the specified resource and returns the result
     * as a JSON object.
     */
    public JsonNode request(String resource, JsonNode data) throws Exception {
        String dataAsString = data.toString();
        return requestImpl(resource, dataAsString, true).responseJson;
    }
}
