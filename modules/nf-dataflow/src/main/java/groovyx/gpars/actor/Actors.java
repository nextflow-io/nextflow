// GPars - Groovy Parallel Systems
//
// Copyright © 2008-2013  The original author or authors
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
// Modified by Seqera Labs for Nextflow (2026): removed the actor factory methods, keeping only the default parallel group, which runs on the dataflow scheduler pool.

package groovyx.gpars.actor;

import groovyx.gpars.group.DefaultPGroup;
import groovyx.gpars.scheduler.Scheduler;

/**
 * Holds the default parallel group shared by actors, agents and dataflow expressions.
 *
 * @author Vaclav Pech, Alex Tkachman
 *         Date: Feb 18, 2009
 */
@SuppressWarnings({"UtilityClass", "AbstractClassWithoutAbstractMethods", "AbstractClassNeverImplemented", "ConstantDeclaredInAbstractClass"})
public abstract class Actors {

    /**
     * The default actor group to share by all actors created through the Actors class.
     */
    public static final DefaultPGroup defaultActorPGroup = new DefaultPGroup(Scheduler.pool());
}
