/*
 * Copyright 2020, Seqera Labs
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
 *
 */

package nextflow.cloud.google.batch.model

import groovy.transform.CompileStatic

/**
 * Model bath lifecycle policy
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class LifecyclePolicy {

    enum Action {
        // Action unspecified.
        ACTION_UNSPECIFIED(0),

        // Action that tasks in the group will be scheduled to re-execute.
        RETRY_TASK(1),

        // Action that tasks in the group will be stopped immediately.
        FAIL_TASK(2);

        private final int value

        Action(int value) { this.value=value }

        static Action valueOf(int value) {
            return forNumber(value);
        }

        /**
         * @param value The numeric wire value of the corresponding enum entry.
         * @return The enum associated with the given numeric wire value.
         */
        static Action forNumber(int value) {
            switch (value) {
                case 0: return ACTION_UNSPECIFIED;
                case 1: return RETRY_TASK;
                case 2: return FAIL_TASK;
                default: return null;
            }
        }

        @Override
        String toString() {
            return this.name()
        }
    }

    static class ActionCondition {
        List<Integer> exitCodes
    }

    // Action to execute when ActionCondition is true.
    Action action

    // Conditions that decide why a task failure is dealt with a specific action.
    ActionCondition action_condition

}
