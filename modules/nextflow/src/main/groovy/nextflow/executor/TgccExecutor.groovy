/*
 * Copyright 2013-2023, Seqera Labs
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
 * Author : Pierre Lindenbaum PhD Institut du Thorax, Nantes, France
 */

package nextflow.executor
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.Files
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
import nextflow.exception.ProcessSubmitException

/**
 This is an Executor to the TGCC  ( https://www-hpc.cea.fr/en/TGCC.html ) 

"TGCC is an infrastructure for scientific high-performance computing and Big Data, able to host petascale supercomputers.
This supercomputing center has been planned to welcome the first French Petascale machine Curie, funded by GENCI for the PRACE Research Infrastructure, and the next generation of Computing Center for Research and Technology (CCRT)."


 
*/
@Slf4j
@CompileStatic
class TgccExecutor extends SlurmExecutor /* TGCC wraps slurm */ {
    static private Pattern TGCC_SUBMIT_REGEX = ~/Submitted Batch Session (\d+)/
    
    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << '-E' << ('\"-D ' + task.workDir + ' --no-requeue\"')
        result << '-r' << getJobNameFor(task).replaceAll("[^A-Za-z_0-9]+","_")
        result << '-o' << quote(task.workDir.resolve(TaskRun.CMD_OUTFILE))
        result << '-e' << quote(task.workDir.resolve(TaskRun.CMD_ERRFILE))
	
	 // the requested partition (a.k.a queue) name
        if( task.config.queue ) {
            result << '-q' << (task.config.queue.toString())
        } else {
	    result << '-q' << 'milan'
        }
	
        if( task.config.project ) {
            result << '-A' << (task.config.project.toString())
        } else {
           throw new UnsupportedOperationException("[TGCC executor] Cannot create a directive without a defined 'project'");
        }

	//  Mounted filesystems
        if( task.config.filesystems ) {
            result << '-m' << (task.config.filesystems.toString())
        } else {
            result << '-m' << 'scratch,store,work,genostore'
        }

        if( task.config.getCpus() > 1 ) {
            result << '-c' << task.config.getCpus().toString()
        }

        if( task.config.time ) {
            	result << '-T' << task.config.getTime().toSeconds().toString()
        	}
        else
        	{
        	result << '-T' << "86400"
        	}

        if( task.config.getMemory() ) {
            result << '-M' << task.config.getMemory().toMega().toString()
        }



        // -- at the end append the command script wrapped file name
        if( task.config.clusterOptions ) {
            result << task.config.clusterOptions.toString() << ''
        }
        
        return result
    }

    @Override
    String getHeaderToken() { '#MSUB' }

    @Override
    protected String getSubmidCmd() {
    	return 'ccc_msub'
    	}
    
    @Override
    protected Pattern  getSubmitRegex() {
    	return 	TGCC_SUBMIT_REGEX;
    }

    @Override
    protected List<String> getKillCommand() { ['ccc_mdel'] }



    @Override
    protected List<String> queueStatusCommand(Object queue) {

        final result = ['sacct','--noheader','-o','jobid%-14s,State%-20s']

        //if( queue )
        //    result << '-p' << queue.toString()

        final user = System.getProperty('user.name')
        if( user )
            result << '-u' << user
        else
            log.debug "[TgccExecutor]:queueStatusCommand Cannot retrieve current user"

        return result
    }

   private QueueStatus decodeQueueStatus(final String s)
   	{
        if(s.equals("COMPLETED")) {
   		return QueueStatus.DONE;
   		}
   	else if(s.equals("RUNNING")) {
   		return QueueStatus.RUNNING;
   		}
   	else if(s.equals("FAILED")) {
   		return QueueStatus.ERROR;
   		}
   	else if(s.equals("NODE_FAIL")) {
   		return QueueStatus.ERROR;
   		}
   	else if(s.equals("PENDING")) {
   		return QueueStatus.PENDING;
   		}
   	else if(s.contains("CANCELLED")) {
   		return QueueStatus.ERROR;
   		}
	else if(s.equals("TIMEOUT")) {
   		return QueueStatus.ERROR;
   		}
	else if(s.equals("OUT_OF_MEMORY")) {
   		return QueueStatus.ERROR;
   		}
	else if(s.equals("SUSPENDED")) {
   		return QueueStatus.HOLD;
   		}
   	else
   		{
   		 log.error "[TGCC Executor] invalid status identifier for Status: `$s` . Interpretted as ERROR. "
   		return QueueStatus.ERROR;
   		}
   	}


    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        def result = [:]

        text.eachLine { String line ->
            def cols = line.split(/\s+/)
            if( cols.size() == 2 ) {
                result.put( cols[0],decodeQueueStatus(cols[1]) )
            	}
            else if(line.contains("CANCELLED"))
            	{
            	result.put( cols[0],  QueueStatus.ERROR )
            	}
            else {
                log.debug "[TGCC Executor] invalid status line: `$line`. Interpreted as ERROR. "
                result.put( cols[0], QueueStatus.ERROR )
            }
        }

        return result
    }


}
