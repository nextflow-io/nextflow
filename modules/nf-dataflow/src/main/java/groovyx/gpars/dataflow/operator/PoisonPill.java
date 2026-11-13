// GPars - Groovy Parallel Systems
//
// Copyright © 2008-2012  The original author or authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Modified by Seqera Labs for Nextflow (2026): a single shared instance; removed the immediate poison pill,
// the countdown hook and the ControlMessage interface.

package groovyx.gpars.dataflow.operator;

/**
 * Represents the poison for dataflow operators.
 * After receiving the poison a dataflow operator will send the poison to all its output channels and terminate.
 *
 * @author Vaclav Pech
 *         Date: Oct 6, 2010
 */
public final class PoisonPill {
    private static final PoisonPill ourInstance = new PoisonPill();

    /**
     * Retrieves the poison pill instance
     * @return The shared singleton poison pill instance
     */
    public static PoisonPill getInstance() {
        return ourInstance;
    }

    private PoisonPill() {
    }
}
