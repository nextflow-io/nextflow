/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.script.control;

/**
 * Extended compilation phases for the Nextflow compiler.
 *
 * @see org.codehaus.groovy.control.Phases
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class Phases {
    public static final int SYNTAX = 1;
    public static final int INCLUDE_RESOLUTION = 2;
    public static final int NAME_RESOLUTION = 3;
    public static final int TYPE_CHECKING = 4;

    public static final int ALL = TYPE_CHECKING;
}
