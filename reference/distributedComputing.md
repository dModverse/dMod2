# Run R Code as Tasks on Other Processes, Machines or a SLURM Cluster

Saves the objects named in `input`, copies them with the code and the
model files to the nodes given by `on` and runs the code there once per
task, in the background. Returns functions to check, wait for, collect
and delete the results. On Windows the local commands run in the sh of
Rtools; remote nodes need a POSIX shell.

## Usage

``` r
distributedComputing(
  expr,
  tasks = NULL,
  on = localNode(),
  jobname = NULL,
  input = ls(parent.frame()),
  build = c("auto", "compile", "link", "copy"),
  libs = NULL,
  seed = NULL,
  wait = FALSE,
  recover = FALSE,
  buildCores = NULL,
  buildBundle = 50
)
```

## Arguments

- expr:

  R code, evaluated once per task.

- tasks:

  A positive whole number of tasks, or a list or vector with one element
  per task, or a data frame with one row per task. `NULL` (default) is
  one task per worker.

- on:

  Where the tasks run:
  [`localNode()`](https://dmodverse.github.io/dMod2/reference/localNode.md)
  (default),
  [`sshNodes()`](https://dmodverse.github.io/dMod2/reference/localNode.md)
  or
  [`slurmCluster()`](https://dmodverse.github.io/dMod2/reference/localNode.md).

- jobname:

  Character, the name of the job and of the job folder. `NULL` (default)
  picks a random name.

- input:

  Character vector, the objects saved and copied to the nodes, looked up
  from the calling environment. Defaults to all objects of the calling
  environment.

- build:

  How the model files reach the nodes: `"auto"`, `"compile"`, `"link"`
  or `"copy"`, see Details.

- libs:

  Optional character vector of library paths put in front of the library
  search path on the nodes, for the tasks and the build. The remote
  shell expands them, so `~` is the remote home.

- seed:

  Integer, task `i` starts with `set.seed(seed + i)`. `NULL` (default)
  draws one.

- wait:

  Logical. `TRUE` waits for the tasks, returns the results of
  [`get()`](https://rdrr.io/r/base/get.html) and deletes the job with
  `purge(local = TRUE)`; a job with failed tasks is kept. Default
  `FALSE`.

- recover:

  Logical. `TRUE` returns the functions of the job `jobname`, submitted
  from this working directory or registered, without submitting again.
  Else, with `on` given, the manifest is taken from the job folder of
  the remote nodes, which recovers a job submitted from another machine;
  it is then saved in `<jobname>_folder` of the working directory and
  registered. A job that also runs on `"localhost"` is recovered only on
  the machine it was submitted from. Default `FALSE`.

- buildCores:

  Number of parallel compiler processes of a build. `NULL` (default)
  uses the cores the node reports, at most 16.

- buildBundle:

  Number of generated sources a build puts into one translation unit,
  only for builds with more files than the shell accepts as arguments.
  `1` compiles one file at a time. Default 50.

## Value

With `wait = TRUE` the results of
[`get()`](https://rdrr.io/r/base/get.html), else an object of class
`dMod2Job`, a list of functions:

- `check()`:

  reports how many tasks are done, running or failed and returns `TRUE`
  once no task runs any more, with a warning when some failed. A task
  fails when the task process ends without a result. A task without a
  result on a node that does not answer, or whose queue does not, counts
  as running.

- `wait(interval = 10)`:

  polls every `interval` seconds until no task runs, reports each change
  as a message and returns `TRUE` when every task is done. While a node
  does not answer, the pause between polls doubles up to 10 minutes.

- `get(combine = TRUE, all = FALSE)`:

  returns the results. With `combine = TRUE` parlists are joined by
  [`c()`](https://rdrr.io/r/base/c.html), fits of
  [`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md) that
  each ran a part of the same starts into one fit, and data frames such
  as parframes by [`rbind()`](https://rdrr.io/r/base/cbind.html), other
  results come as a list; tasks without a result are left out with a
  warning. `combine = FALSE` returns one entry per task, named like
  `tasks`, `NULL` for a task not done. Results of remote nodes are
  copied into `<jobname>_folder/<k>`, with `all = TRUE` together with
  every other file of the job folder.

- `log(i, lines = 20)`:

  prints the last lines of the output of task `i`.

- `purge(local = FALSE)`:

  deletes the job folders on the nodes, with `local = TRUE` also the
  local `<jobname>_folder` and the registry entry.

- `terminate()`:

  stops the tasks.

## Details

The start is reported as a message, on a SLURM cluster with the job ID.
A build writes the output to `build.log` in the job folder; the error of
a failed build shows the end of that log.

The value of `expr` is the result of a task. In the code, `task` is the
task's element of `tasks` (`NULL` for a count), `taskId` the position of
the task, `cores` the cores the task may use and `jobname` the job's
name. Parallel code such as
[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md) or
[profile()](https://dmodverse.github.io/dMod2/reference/profile.md)
takes `cores = cores`. A task that raises an error returns the error as
a `try-error`. Task `i` starts with the seed `seed + i`.

The tasks are dealt in turn to the workers of `on`, the nodes or the
array tasks of a SLURM job, and a worker runs the tasks dealt to it one
after the other, see
[`localNode()`](https://dmodverse.github.io/dMod2/reference/localNode.md).
A job is kept in `<jobname>_folder`, in the working directory for
`"localhost"` and in the home directory of a remote node. The job
manifest `job.rds` lies in the local `<jobname>_folder` and in that of
every remote node. A job of the same name that still runs there is an
error, also when it was submitted from another machine. Every job is
registered in `tools::R_user_dir("dMod2", "data")` until
`purge(local = TRUE)`, and
[`dMod2Jobs()`](https://dmodverse.github.io/dMod2/reference/dMod2Jobs.md)
lists them.

The model files are those of the
[obsfn](https://dmodverse.github.io/dMod2/reference/obsfn.md),
[parfn](https://dmodverse.github.io/dMod2/reference/parfn.md),
[prdfn](https://dmodverse.github.io/dMod2/reference/prdfn.md) and
objective objects in `input`; objects without them take every file of
`getOption("dMod2.outdir")`, else of the working directory.

- `build = "compile"`: the C/C++ sources are copied and compiled once
  per node into one shared object, against the cppDE installed there.

- `build = "link"`: the object files (`.o`) are copied and linked once
  per node; needs a compiler there compatible with the local one.

- `build = "copy"`: the shared objects are copied.

- `build = "auto"` (default): `"copy"` on the local machine, `"compile"`
  elsewhere.

After a build the
[obsfn](https://dmodverse.github.io/dMod2/reference/obsfn.md),
[parfn](https://dmodverse.github.io/dMod2/reference/parfn.md) and
[prdfn](https://dmodverse.github.io/dMod2/reference/prdfn.md) objects
and the objectives use the new shared object.

## See also

[`localNode()`](https://dmodverse.github.io/dMod2/reference/localNode.md),
[`slurmTemplate()`](https://dmodverse.github.io/dMod2/reference/slurmTemplate.md),
[`dMod2Jobs()`](https://dmodverse.github.io/dMod2/reference/dMod2Jobs.md),
[`detectFreeCores()`](https://dmodverse.github.io/dMod2/reference/detectFreeCores.md);
the guide chapter [Remote
computing](https://dmodverse.github.io/dMod2/articles/guide-remote-computing.html).

## Examples

``` r
# \donttest{
# Two tasks on the local machine, in a temporary working directory; each task
# starts an R process
if (.Platform$OS.type == "unix") {
  old <- setwd(tempdir())
  M <- matrix(rnorm(1e2), 10, 10)
  job <- distributedComputing(det(M * task), tasks = list(a = 1, b = 2),
                              input = "M", jobname = "example")
  job$wait(interval = 0.1)
  print(job$get())
  job$purge(local = TRUE)
  setwd(old)
}
#> example: 2 tasks started on localhost (3 cores)
#> example: 0 of 2 tasks done, 2 running
#> example: 1 of 2 tasks done, 1 running
#> example: 2 of 2 tasks done
#> $a
#> [1] -403.4075
#> 
#> $b
#> [1] -413089.3
#> 
# }

if (FALSE) { # \dontrun{
# For an objective obj and a center of the starts: multistarts on two nodes
# reached by ssh, one task per node with ten fits per core
nodes <- sshNodes(c("server1", "server2"), cores = 16)
jobFitting <- distributedComputing(
  mstrust(obj, center, fits = 10 * cores, cores = cores, sd = 4),
  on = nodes, jobname = "fits"
)
jobFitting$wait()
bestfit <- as.parvec(jobFitting$get())
jobFitting$purge(local = TRUE)

# Profiles of the best fit: one task per group of four parameters and side, on
# the free cores of both nodes. get() joins the sides of every parameter.
groups <- split(names(bestfit), ceiling(seq_along(bestfit) / 4))
jobProfiling <- distributedComputing(
  profile(obj, bestfit, whichPar = groups[[task$group]], side = task$side,
          limits = c(-5, 5), cores = cores),
  tasks = expand.grid(group = seq_along(groups), side = c("left", "right"),
                      stringsAsFactors = FALSE),
  on = sshNodes(c("server1", "server2")), jobname = "profiles"
)
jobProfiling$wait()
profiles <- jobProfiling$get()
confint(profiles)

# 1000 fits on a SLURM cluster: 100 array tasks of 10 fits on 10 cores each,
# through the job script template of the login node "mycluster", set up once
# from a job script of the cluster with the placeholders of ?slurmTemplate.
slurmTemplate("mycluster", from = "~/mycluster.sh")
cluster <- slurmCluster("mycluster", cores = 10, walltime = "02:00:00")
job <- distributedComputing(
  mstrust(obj, center, fits = cores, cores = cores, sd = 4),
  tasks = 100, on = cluster, jobname = "slurmfits"
)

# The ten starts of one mixed-effects multistart as ten array tasks of four cores:
# every task draws the same starts and fits a separate part, get() joins them into one fit
jobEm <- distributedComputing({
  set.seed(1)
  msem(obj, random, init, fits = 10, which = taskId,
       cores = c(fits = 1, conditions = cores))
}, tasks = 10, on = slurmCluster("mycluster", cores = 4, walltime = "06:00:00"),
jobname = "msemfits")

# The jobs submitted and not purged, and the functions of one of them
dMod2Jobs()
job <- distributedComputing(NULL, jobname = "slurmfits", recover = TRUE)
} # }
```
