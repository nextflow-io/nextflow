(metrics-page)=

# Understanding task resource metrics

This tutorial explains how the resource usage metrics from the {ref}`Execution report <execution-report>` are computed. See {ref}`Execution report <execution-report>` for more information about how to generate an execution report.

## CPU Usage

The CPU Usage plot reports how CPU resources are used by each process.

```{image} _static/report-resource-cpu.png
```

The **Raw Usage** tab is expected to show 100% core usage if a process performs one task of pure computation. If the task is distributed over, 2, 3, or 4 CPUs, the raw usage will be 200%, 300%, or 400%, respectively. The **% Allocated** tab rescales the raw usage value relative to the number of CPUs that are set with the `cpus` directive. If the `cpus` directive is not set, the number of CPUs is set to `1` and the **% Allocated** tab will show the same values the **Raw Usage** tab.

For example, using the program [stress](https://people.seas.harvard.edu/~apw/stress/), the following script would report 100% CPU usage in the **Raw Usage** tab and 50% CPU usage in the **% Allocated** tab as the process requested double the number of CPUs that are required:

```nextflow
process cpuUsageEx1 {
  cpus 2

  script:
  """
  stress -c 1 -t 10 # compute square-root of random numbers during 10s using 1 CPU
  """
}

workflow{
    cpuUsageEx1()
}
```

The CPU usage decreases if the process spends some time performing pure computation and some time waiting the CPU. For example, using the program [stress](https://people.seas.harvard.edu/~apw/stress/) and `sleep`, the following script would report 75% CPU usage in the **Raw Usage** tab:

```nextflow
process cpuUsageEx2 {
  cpus 1

  script:
  """
  stress -c 1 -t 10 # compute the square-root of random numbers during 10s using 1 CPU
  stress -c 1 -t 5 # compute the square-root of random numbers during 5s using 1 CPU
  sleep 5 # use no CPU during 5s
  """
}

workflow{
    cpuUsageEx2()
}
```

In the above example, the CPU usage is a weighted average that accounts for the percentage of the CPU used and duration of each individual program over the job duration (elapsed real time, real time, or wall time):

$$
\frac{ 100\% \times 10s + 100\% \times 5s + 0\% \times 5s }{10s+5s+5s} = 75\%
$$

The CPU usage increases if a single step is forked on multiple CPUs:

```nextflow
process cpuUsageEx3 {
  cpus 2

  script:
  """
  stress -c 2 -t 10 # compute square-root of random numbers during 10s using 2 CPUs
  sleep 10 # use no CPU during 10s
  """
}

workflow{
    cpuUsageEx3()
}
```

In the above example, the **Raw Usage** tab would report 100%:

$$
\frac{ 200\% \times 10s + 0\% \times 10s }{10s+10s} = 100\%
$$

However, the **% Allocated** tab would report 50%. It would not be relevant to change the `cpus` directive from `2` to `1` as the process uses 2 CPUs at it peak load.

:::{tip}
The [stress](https://people.seas.harvard.edu/~apw/stress/) program can be installed with `sudo apt-get install stress` or `sudo yum install stress` depending on your Linux distribution.
:::

## Memory Usage

The Memory Usage plot reports how memory was used by each process. It has three tabs, **Physical (RAM)**, **Virtual (RAM + Disk swap)**, and **% RAM Allocated**, showing the usage of the physical memory (RAM), the virtual memory (`vmem`), and the percentage of RAM used by the process relative to the memory that the `memory` directive set, respectively.

The peak usage during the execution of the process is reported for both physical and virtual memories. The total amount of memory used by a process is the `virtual memory (vmem)`. The `vmem` contains all memory areas, including in the physical memory (RAM), in the swap space, on the disk, or shared with other processes. The `resident set size (RSS)` is the amount of `physical memory (RAM)` held by a process.

The relationship is:

$$
vmem \geq RSS + Swap
$$

The behavior of the memory usage plot can be examined using two programs written in C. The first program allocates a variable of 1 GiB:

```{code-block} c
:emphasize-lines: 31,43
:linenos: true

#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>

/* Get vmem and rss usage from /proc/<pid>/statm */
static int mem_used(pid_t pid, unsigned long* vmem, unsigned long* rss) {
    FILE* file;
    char path[40];
    unsigned int page_size;

    snprintf(path, 40, "/proc/%ld/statm", (long) pid);
    file = fopen(path, "r");
    // vmem and rss are the first values in the file
    fscanf(file, "%lu %lu", vmem, rss);
    // values in statm are in pages so to get bytes we need to know page size
    page_size = (unsigned) getpagesize();
    *vmem = *vmem * page_size;
    *rss = *rss * page_size;

    fclose(file);
    return 0;
}

int main(int argc, char **argv) {
    unsigned char *address;
    char input;
    size_t size = 1024*1024*1024;  // 1 GiB
    unsigned long i;
    unsigned long vmem = 0;
    unsigned long rss = 0;
    pid_t pid;

    pid = getpid();
    printf("Pid: %ld\n", (long) pid);

    mem_used(pid, &vmem, &rss);
    printf("VMEM: %lu RSS: %lu\n", vmem, rss);

    address = malloc(size);
    printf("Allocated %d Bytes of memory\n", (int) size);

    mem_used(pid, &vmem, &rss);
    printf("VMEM: %lu RSS: %lu\n", vmem, rss);

    // Leave time for nextflow to get information
    sleep(15);

    free(address);
    return 0;
}
```

The second program allocates a variable of 1 GiB and fills it with data:

```{code-block} c
:emphasize-lines: 31,43,49-53
:linenos: true

#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>

/* Get vmem and rss usage from /proc/<pid>/statm */
static int mem_used(pid_t pid, unsigned long* vmem, unsigned long* rss) {
    FILE* file;
    char path[40];
    unsigned int page_size;

    snprintf(path, 40, "/proc/%ld/statm", (long) pid);
    file = fopen(path, "r");
    // vmem and rss are the first values in the file
    fscanf(file, "%lu %lu", vmem, rss);
    // values in statm are in pages so to get bytes we need to know page size
    page_size = (unsigned) getpagesize();
    *vmem = *vmem * page_size;
    *rss = *rss * page_size;

    fclose(file);
    return 0;
}

int main(int argc, char **argv) {
    unsigned char *address;
    char input;
    size_t size = 1024*1024*1024;  // 1 GiB
    unsigned long i;
    unsigned long vmem = 0;
    unsigned long rss = 0;
    pid_t pid;

    pid = getpid();
    printf("Pid: %ld\n", (long) pid);

    mem_used(pid, &vmem, &rss);
    printf("VMEM: %lu RSS: %lu\n", vmem, rss);

    address = malloc(size);
    printf("Allocated %d Bytes of memory\n", (int) size);

    mem_used(pid, &vmem, &rss);
    printf("VMEM: %lu RSS: %lu\n", vmem, rss);

    printf("Filling memory with data...");
    fflush(stdout);
    for (i = 0; i < size; i++) {
        *(address + i) = 123;
    }

    mem_used(pid, &vmem, &rss);
    printf("\nVMEM: %lu RSS: %lu\n", vmem, rss);

    // Leave time for nextflow to get information
    sleep(15);

    free(address);
    return 0;
}
```

The first and second programs are executed as foo` and `bar` in the following script:

```nextflow
process foo {
    memory '1.5 GB'

    script:
    """
    memory_vmem_1GiB_ram_0Gib
    """
}

process bar {
    memory '1.5 GB'

    script:
    """
    memory_vmem_1GiB_ram_1Gib
    """
}

workflow{
    foo() // Allocates a variable of 1 GiB
    bar() // Allocates a variable of 1 GiB and fills it with data
}
```

The **Virtual (RAM + Disk swap)** tab shows that both `foo` and `bar` use the same amount of virtual memory (~1 GiB):

```{image} _static/report-resource-memory-vmem.png
```

However, the **Physical (RAM)** tab shows that the `bar` uses ~1 GiB of RAM while `foo` uses ~0 GiB:

```{image} _static/report-resource-memory-ram.png
```

The **% RAM Allocated** tab shows that 0% of the resource set in the `memory` directive was used for `foo` and 67% of the resource are used for `bar`:

```{image} _static/report-resource-memory-pctram.png
```

:::{warning}
Memory and storage metrics are reported in bytes. For example, 1 KB = $1024$ bytes, 1 MB = $1024^2$ bytes, and 1 GB = $1024^3$ bytes.
:::

## Job Duration

The Job Duration plot reports how long each process took to run. It has two tabs. The **Raw Usage** tab shows the job duration (a.k.a. elapsed real time, real time or wall time ) and the **% Allocated** tab shows the time that was requested relative to what was requested using the `time` directive.

```{image} _static/report-resource-job-duration.png
```

## I/O Usage

The I/O Usage plot shows how much data was read and written by a processes.

The amount of data that was read by a process (`rchar` in the trace file) is the number of bytes the process read, using any read-like system calls. For example, from files, pipes, and tty. The amount of data that was written by a process (`wchar` in the trace file) is the number of bytes the process wrote, using any write-like system call.

The I/O plot **Read** tab shows how much data was read and the **Write** tab shows how much data was written by each process. For example, the following script reads and writes different data volumes:

```nextflow
process io_read_write_1G {
    script:
    """
    dd if=/dev/zero of=/dev/null bs=1G count=1
    """
}

process io_read_write_256M {
    script:
    """
    dd if=/dev/zero of=/dev/null bs=256M count=1
    """
}

workflow{
    io_read_write_1G() // Read and write 1 GiB
    io_read_write_256M() // Read and write 256 Mb
}
```

The **Read** tab shows that ~1 Gib and ~256 Mb are read:

```{image} _static/report-resource-io-read.png
```

The **Write** tab shows that ~1 Gib and ~256 Mb are written:

```{image} _static/report-resource-io-write.png
```
