# Distributed computing

`distributedComputing()` runs R code as tasks elsewhere: as background processes on
the local machine, on machines reached by ssh, or as an array job of a SLURM cluster. It
saves the objects the code needs, copies them with the code to every machine, builds
the models there if asked, and returns functions to check, collect and delete the
results. The tasks run detached from the R session that submitted them.
Two demos use it. [Bile acid transport](BA_transport.md) runs its multistarts as tasks
on the local machine.
[IFN-alpha sensitization](KokMSB2020.md) ran its 200 fits as a SLURM array job on the
bwForCluster Helix; without an account there the demo reads the fits it ships.

## Machines, cores and tasks

`machine` names the machines, each once: as a character vector, or as a vector of
cores named by machine.

```r
machine = c("server1", "server2")           # every machine with its free cores
machine = c(server1 = 20, server2 = 10)     # with these cores
```

`"localhost"` is the local machine, used without ssh, with the job folder in the working
directory. Any other name is an ssh destination such as `"host"` or `"user@host"`,
reached with key authentication, with the job folder in its home directory. Before the
start `detectFreeCores()` looks at every machine; a machine gets at most its free cores,
with a warning when it was given more.

The code runs once per task. By default there is one task per machine. `tasks` is
otherwise a number of tasks, or a list or data frame with one element or row per task.
The tasks are dealt to the machines in turn, and a machine runs its tasks one after the
other, each with all its cores. The code reads its element or row as `task`, its
position as `taskId` and the cores it may use as `cores`, which parallel code such as
`mstrust()` or `profile()` takes as it is. Every task starts with the seed
`seed + taskId`, and an error in a task becomes its result as a `try-error`.

For an objective `obj` and a centre `center` of the starts, multistarts and then the
profiles of the best fit:

```r
# one task per machine, ten multistarts per core: 300 fits
jobFitting <- distributedComputing(
  mstrust(obj, center, fits = 10 * cores, cores = cores),
  machine = c(server1 = 20, server2 = 10), jobname = "fits"
)
while (!jobFitting$check()) Sys.sleep(10)
parlists <- Filter(function(f) inherits(f, "parlist"), jobFitting$get())
fits     <- as.parframe(do.call(c, parlists))
bestfit  <- as.parvec(fits)

# profiles of the best fit, four parameters per task, on the free cores of both machines
jobProfiling <- distributedComputing(
  profile(obj, bestfit, whichPar = task, cores = cores),
  tasks = split(names(bestfit), ceiling(seq_along(bestfit) / 4)),
  machine = c("server1", "server2"), jobname = "profiles"
)
while (!jobProfiling$check()) Sys.sleep(10)
profiles <- do.call(rbind, jobProfiling$get())
plotProfile(profiles)
```

The objects of `input`, by default every object of the calling environment, travel
with the code; the code attaches dMod2 and nothing else.

## Model files and builds

The model files come from `getOption("dMod2.outdir")`, else from the working directory. Without `compile` and `link` the shared objects are copied, which needs a machine
compatible with the local one. With `compile = TRUE` the C and C++ sources are copied and
compiled once per machine into one shared object, against the cppDE installed there;
`link = TRUE` copies and links the object files. `libs` puts library paths first on
the remote side, for the build and the tasks.

## Following a job

```r
jobFitting$check()               # tasks done, running and failed; TRUE once none runs
results <- jobFitting$get()      # one result per task, NULL for a task not done
jobFitting$get(all = TRUE)       # also every other file of the job folders
jobFitting$terminate()           # stops the tasks
jobFitting$purge(local = TRUE)   # deletes the job folders, and the local one

# the functions of the job "fits" submitted from this working directory, in a new session
jobFitting <- distributedComputing(NULL, jobname = "fits", recover = TRUE)
```

`wait = TRUE` waits for the tasks, returns the results and deletes the job.

## A SLURM cluster

With `slurm = TRUE` the tasks form an array job, submitted on the login node named by
`machine` through the job script template set up for it, and `machine` names the cores
of each array task, such as `c(mycluster = 20)`; the code reads them as `cores`.

The template is set up once per login node with `slurmTemplate()` and kept in the
configuration directory of dMod2, `tools::R_user_dir("dMod2", "config")`, across
sessions. It is a job script with the fields `{{ jobname }}`, `{{ ntasks }}`,
`{{ cores }}`, `{{ memPerCore }}` and `{{ walltime }}` and one line `{{ run }}`, which
runs the task of `SLURM_ARRAY_TASK_ID`. The lines before `{{ run }}` that are not
comments set up the environment, for example with `module load`, and the build on the
login node runs in it too. dMod2 ships a template for the bwForCluster Helix.

```r
slurmTemplate("mycluster", from = "~/mycluster.sh")   # once, opens it for editing
slurmTemplate("helix")                                 # the template for Helix
slurmTemplate()                                        # the templates set up

job <- distributedComputing(mstrust(obj, center, fits = 100, cores = cores),
                            tasks = 4, machine = c(mycluster = 20), slurm = TRUE,
                            walltime = "02:00:00", compile = TRUE)
```

`jobscript` passes another template for one job.

## Access

A remote machine needs ssh access with a key, without a password prompt, and dMod2 and
cppDE installed for the R that runs the tasks there; `libs` points to that library when
it is not on the default path. A SLURM cluster needs an account on its login node and a
job script adapted to it.

## Requirements

A Unix-like system on every side, ssh with key authentication to every remote machine,
and dMod2 with cppDE installed for the R that runs the tasks there. `zstd`, where
present locally, compresses the transferred workspace and is then needed on the remote
side too. `walltime` ends a task without SLURM through `timeout` where the machine
has it.
