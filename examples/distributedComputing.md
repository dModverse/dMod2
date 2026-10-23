# Distributed computing

`distributedComputing()` runs R code as tasks elsewhere: as background processes on
the local machine, on machines reached by ssh, or as an array job of a SLURM cluster. It
saves the objects the code needs, copies them with the code and the model files,
builds the models there if needed, and returns functions to follow, collect and delete
the results. The tasks run detached from the R session that submitted them.
Two demos use it. [Bile acid transport](BA_transport.md) runs its multistarts as tasks
on the local machine.
[IFN-alpha sensitization](KokMSB2020.md) ran its 200 fits as a SLURM array job on the
bwForCluster Helix; without an account there the demo reads the fits it ships.

## Nodes, tasks and cores

`on` says where the tasks run, and every task gets `cores` cores:

```r
on = localNode()                                         # the local machine, its free cores
on = localNode(cores = 8)
on = sshNodes(c("server1", "server2"))                   # every node with its free cores
on = sshNodes(c("server1", "server2"), cores = 10)
on = sshNodes(cores = c(server1 = 24, server2 = 20))
on = slurmCluster("mycluster", cores = 8, walltime = "12:00:00")
on = slurmCluster("mycluster", cores = 8, array = 32, walltime = "12:00:00")
```

A node is an ssh destination such as `"host"` or `"user@host"`, reached with key
authentication, or `"localhost"`, the local machine without ssh. Before the start
`detectFreeCores()` looks at every node; a node gets at most its free cores, with a
warning when it was given more. A SLURM cluster is reached through its login node, and
`array` is the size of its array job, by default one array task per task.

The tasks run on workers: a node is one worker, an array task another. `tasks` is a
number of tasks, or a list or data frame with one element or row per task; by default
there is one task per worker. The tasks are dealt to the workers in turn, and a worker
runs its tasks one after the other. The code reads its element or row as `task`, its
position as `taskId` and its cores as `cores`, which parallel code such as `mstrust()`
or `profile()` takes as it is. Every task starts with the seed `seed + taskId`, and an
error in a task becomes its result as a `try-error`.

For an objective `obj` and a centre `center` of the starts, multistarts and then the
profiles of the best fit:

```r
# one task per node, ten multistarts per core: 300 fits
jobFitting <- distributedComputing(
  mstrust(obj, center, fits = 10 * cores, cores = cores),
  on = sshNodes(cores = c(server1 = 20, server2 = 10)), jobname = "fits"
)
jobFitting$wait()
fits    <- as.parframe(jobFitting$get())
bestfit <- as.parvec(fits)

# profiles of the best fit, one task per group of four parameters and side
groups <- split(names(bestfit), ceiling(seq_along(bestfit) / 4))
jobProfiling <- distributedComputing(
  profile(obj, bestfit, whichPar = groups[[task$group]], side = task$side, cores = cores),
  tasks = expand.grid(group = seq_along(groups), side = c("left", "right"),
                      stringsAsFactors = FALSE),
  on = sshNodes(c("server1", "server2")), jobname = "profiles"
)
jobProfiling$wait()
profiles <- jobProfiling$get()
plotProfile(profiles)
```

`get()` joins the results of the tasks: parlists with `c()`, data frames and parameter
frames with `rbind()`, which also joins the two sides of a profile into one. A task
without a result is left out with a warning; `get(combine = FALSE)` returns one entry per
task.

The objects of `input`, by default every object of the calling environment, travel
with the code; the code attaches dMod2 and nothing else.

## Model files and builds

The model files are those of the model objects in `input`: their C and C++ sources,
object files or shared objects. Objects without model files take every such file of
`getOption("dMod2.outdir")`, else of the working directory. `build` says how they reach
the nodes:

- `"auto"` (default): `"copy"` on the local machine, `"compile"` elsewhere.
- `"compile"`: the sources are compiled once per node into one shared object, against
  the cppDE installed there.
- `"link"`: the object files are linked once per node, which needs a compiler there
  compatible with the local one.
- `"copy"`: the shared objects are copied, which needs a node compatible with the local
  machine.

`libs` puts library paths first on the remote side, for the build and the tasks.

## Following a job

```r
jobFitting$check()               # tasks done, running and failed; TRUE once none runs
jobFitting$wait()                # waits until none runs, reporting each change
results <- jobFitting$get()      # the results joined
jobFitting$get(all = TRUE)       # also every other file of the job folders
jobFitting$log(3)                # the last lines of the output of task 3
jobFitting$terminate()           # stops the tasks
jobFitting$purge(local = TRUE)   # deletes the job folders, and the local one

dMod2Jobs()                      # the jobs submitted and not purged
# the functions of the job "fits", in a new session
jobFitting <- distributedComputing(NULL, jobname = "fits", recover = TRUE)
```

`wait = TRUE` waits for the tasks, returns the results and deletes the job. A job keeps
its name while it runs: submitting another job of that name where it runs is an error.
Every job is registered in the data directory of dMod2,
`tools::R_user_dir("dMod2", "data")`, until `purge(local = TRUE)`, so `recover` finds it
from any working directory.

## A SLURM cluster

On a `slurmCluster()` every worker is an array task with `cores` cores, submitted on the
login node through the job script template set up for it. `partition`, `account`,
`maxRunning` (the most array tasks at once) and `sbatch` (further options) go to
`sbatch` and take precedence over the `#SBATCH` lines of the template.

The template is set up once per login node with `slurmTemplate()` and kept in the
configuration directory of dMod2, `tools::R_user_dir("dMod2", "config")`, across
sessions. It is a job script with the fields `{{ jobname }}`, `{{ ntasks }}` (the number
of array tasks), `{{ cores }}`, `{{ memPerCore }}` and `{{ walltime }}` and one line
`{{ run }}`, which runs the tasks of `SLURM_ARRAY_TASK_ID`. The lines before `{{ run }}`
that are not comments set up the environment, for example with `module load`, and the
build on the login node runs in it too. dMod2 ships a template for the bwForCluster
Helix.

```r
slurmTemplate("mycluster", from = "~/mycluster.sh")   # once, opens it for editing
slurmTemplate("helix")                                 # the template for Helix
slurmTemplate()                                        # the templates set up

# 1000 fits: 100 array tasks of 10 fits on 10 cores each
job <- distributedComputing(mstrust(obj, center, fits = cores, cores = cores),
                            tasks = 100, jobname = "slurmfits",
                            on = slurmCluster("mycluster", cores = 10, walltime = "02:00:00"))
```

`slurmCluster(jobscript = )` passes another template for one job.

## Requirements

A Unix-like system on every side, ssh with key authentication to every remote node, and
dMod2 with cppDE installed for the R that runs the tasks there; `libs` points to that
library when it is not on the default path. A SLURM cluster needs an account on its login
node and a job script template adapted to it. `zstd`, where present locally, compresses
the transferred workspace and is then needed on the remote side too.
