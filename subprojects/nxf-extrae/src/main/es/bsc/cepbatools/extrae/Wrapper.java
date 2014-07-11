///****************************************************************************\
// *                        ANALYSIS PERFORMANCE TOOLS                         *
// *                                   Extrae                                  *
// *              Instrumentation package for parallel applications            *
// *****************************************************************************
// *     ___     This library is free software; you can redistribute it and/or *
// *    /  __         modify it under the terms of the GNU LGPL as published   *
// *   /  /  _____    by the Free Software Foundation; either version 2.1      *
// *  /  /  /     \   of the License, or (at your option) any later version.   *
// * (  (  ( B S C )                                                           *
// *  \  \  \_____/   This library is distributed in hope that it will be      *
// *   \  \__         useful but WITHOUT ANY WARRANTY; without even the        *
// *    \___          implied warranty of MERCHANTABILITY or FITNESS FOR A     *
// *                  PARTICULAR PURPOSE. See the GNU LGPL for more details.   *
// *                                                                           *
// * You should have received a copy of the GNU Lesser General Public License  *
// * along with this library; if not, write to the Free Software Foundation,   *
// * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA          *
// * The GNU LEsser General Public License is contained in the file COPYING.   *
// *                                 ---------                                 *
// *   Barcelona Supercomputing Center - Centro Nacional de Supercomputacion   *
//\****************************************************************************/
//
///* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*\
// | @file: $HeadURL: https://svn.bsc.es/repos/ptools/extrae/branches/2.4/src/tracer/wrappers/API/wrapper.c $
// | @last_commit: $Date: 2013-11-26 10:30:20 +0100 (Tue, 26 Nov 2013) $
// | @version:     $Revision: 2336 $
//\* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
package es.bsc.cepbatools.extrae;

public class Wrapper
{

    static
    {
        System.loadLibrary("javaseqtrace");
    }

    public static int EXTRAE_DISABLE_ALL_OPTIONS = 0;
    public static int EXTRAE_CALLER_OPTION       = 1 << 0;
    public static int EXTRAE_HWC_OPTION          = 1 << 1;
    public static int EXTRAE_MPI_HWC_OPTION      = 1 << 2;
    public static int EXTRAE_MPI_OPTION          = 1 << 3;
    public static int EXTRAE_OMP_OPTION          = 1 << 4;
    public static int EXTRAE_OMP_HWC_OPTION      = 1 << 5;
    public static int EXTRAE_UF_HWC_OPTION       = 1 << 6;
    public static int EXTRAE_PTHREAD_OPTION      = 1 << 7;
    public static int EXTRAE_PTHREAD_HWC_OPTION  = 1 << 8;
    public static int EXTRAE_SAMPLING_OPTION     = 1 << 9;
    public static int EXTRAE_ENABLE_ALL_OPTIONS  =
            EXTRAE_CALLER_OPTION |
                    EXTRAE_HWC_OPTION |
                    EXTRAE_MPI_HWC_OPTION |
                    EXTRAE_MPI_OPTION |
                    EXTRAE_OMP_OPTION |
                    EXTRAE_OMP_HWC_OPTION |
                    EXTRAE_UF_HWC_OPTION |
                    EXTRAE_PTHREAD_OPTION |
                    EXTRAE_PTHREAD_HWC_OPTION |
                    EXTRAE_SAMPLING_OPTION;

    public static native void Init();
    public static native void Fini();
    public static native void Event(int id, long val);
    public static native void nEvent(int types[], long values[]);
    public static native void defineEventType(int type, String description,
                                              long nValues, long[] values, String[] descriptionValues);
    public static native void SetOptions(int options);
    public static native void Shutdown();
    public static native void Restart();

    public static native void Comm(boolean send, int tag, int size,
                                   int partner, long id);
    public static native int GetPID();
    public static native void SetTaskID(int id);
    public static native void SetNumTasks(int num);
    public static native int GetTaskID();
    public static native int GetNumTasks();
    public static native void SetThreadID(int id);
    public static native void SetNumThreads(int num);
    public static native int GetThreadID();
    public static native int GetNumThreads();

    public static native void resumeVirtualThread(long vthread);
    public static native void suspendVirtualThread();
}

