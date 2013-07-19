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

import com.fasterxml.jackson.core.JsonFactory;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.MappingJsonFactory;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.BooleanNode;
import com.fasterxml.jackson.databind.node.DoubleNode;
import com.fasterxml.jackson.databind.node.IntNode;
import com.fasterxml.jackson.databind.node.LongNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import com.fasterxml.jackson.databind.node.TextNode;
import java.io.IOException;

/**
 * Utility class for working with JSON objects.
 */
public class DXJSON {

    private static final JsonFactory dxJsonFactory = new MappingJsonFactory();
    private static ObjectMapper mapper = new ObjectMapper();

    // Utility class should not be instantiated
    private DXJSON() {
    }

    /**
     * Parses the specified string into a JSON object.
     */
    public static JsonNode parseJson(String stringified) throws IOException {
        return dxJsonFactory.createJsonParser(stringified).readValueAsTree();
    }

    // TODO: helpers for making arrays

    /**
     * Builder class that generates a JSON object (hash).
     *
     * Example:
     *
     * <pre>
     * {@code
     * ObjectNode o = DXJSON.getObjectBuilder()
     *                      .put("key1", "a-string")
     *                      .put("key2", 12321)
     *                      .build()}</pre>
     *
     * when serialized, produces the JSON object <tt>{"key1": "a-string", "key2": 12321}</tt>.
     */
    public static class ObjectBuilder {
        private final boolean isEmpty;
        private final ObjectBuilder next;
        private final String key;
        private final JsonNode value;

        private ObjectBuilder(boolean isEmpty, ObjectBuilder next, String key, JsonNode value) {
            this.isEmpty = isEmpty;
            this.next = next;
            this.key = key;
            this.value = value;
        }

        /**
         * Initializes an ObjectBuilder which will generate an empty object.
         */
        public ObjectBuilder() {
            this(true, null, null, null);
        }

        /**
         * Adds a key-value pair with an arbitrary JsonNode value.
         */
        public ObjectBuilder put(String key, JsonNode value) {
            return new ObjectBuilder(false, this, key, value);
        }

        // TODO: allow easy creation of nulls

        /**
         * Adds a key-value pair with a string value and returns the resulting
         * ObjectBuilder.
         */
        public ObjectBuilder put(String key, String value) {
            return put(key, new TextNode(value));
        }

        /**
         * Adds a key-value pair with a numeric value and returns the resulting
         * ObjectBuilder.
         */
        public ObjectBuilder put(String key, int value) {
            return put(key, new IntNode(value));
        }

        /**
         * Adds a key-value pair with a numeric value and returns the resulting
         * ObjectBuilder.
         */
        public ObjectBuilder put(String key, long value) {
            return put(key, new LongNode(value));
        }

        /**
         * Adds a key-value pair with a numeric value and returns the resulting
         * ObjectBuilder.
         */
        public ObjectBuilder put(String key, double value) {
            return put(key, new DoubleNode(value));
        }

        /**
         * Adds a key-value pair with a boolean value and returns the resulting
         * ObjectBuilder.
         */
        public ObjectBuilder put(String key, boolean value) {
            return put(key, value ? BooleanNode.TRUE : BooleanNode.FALSE);
        }

        /**
         * Generates a JSON object.
         */
        public ObjectNode build() {
            ObjectNode output = mapper.createObjectNode();
            ObjectBuilder nextBuilder = this;
            while (!nextBuilder.isEmpty) {
                output.put(nextBuilder.key, nextBuilder.value);
                nextBuilder = nextBuilder.next;
            }
            return output;
        }
    }

    public static ObjectBuilder getObjectBuilder() {
        return new ObjectBuilder();
    }

}
