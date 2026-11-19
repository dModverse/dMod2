# 11. Remote computing

dMod2 parallelizes along two axes: conditions of one objective
evaluation in threads, independent fits and profiles in processes.
[`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md)
runs any R code as tasks in processes detached from the session,
locally, over ssh or as a SLURM array job. The local job below runs at
page build; commands for other machines are shown, not run.

## Threads and processes

`options(dMod2.cores = n)`, or `obj(pars, cores = n)` for one call,
solves the conditions of an evaluation in `n` threads of an OpenMP team
inside cppDE, constant in size within a session. The guide model has two
conditions, so this axis uses at most two threads; many conditions, or
[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) with many
groups, gain in proportion. A single condition never uses it.
Predictions do not pass through forked processes: sending back the
arrays of states and sensitivities costs more than computing them.

[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md),
[`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md) and
[`profile()`](https://dmodverse.github.io/dMod2/reference/profile.md)
run independent fits as processes through `cores`. With one number, fits
are forks of the session, each solving conditions in one thread: a
forked child never starts the thread team, and BLAS and OpenMP stay at
one thread across the fork. `cores = c(fits = 4, conditions = 2)` runs
four fits with two threads each in fresh R processes (a PSOCK cluster),
as always on Windows. Each process holds a copy of the predictions, so
memory grows with the process count; keep the product of both below the
core count.

| where | argument | axis |
|----|----|----|
| one evaluation | `options(dMod2.cores = )`, `obj(pars, cores = )` | conditions, threads |
| [`em()`](https://dmodverse.github.io/dMod2/reference/em.md) | `cores` | groups and conditions, threads |
| [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md), [`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md) | `cores`, `c(fits = , conditions = )` | fits, processes |
| [`profile()`](https://dmodverse.github.io/dMod2/reference/profile.md) | `cores`, `c(pars = , conditions = )` | parameters, processes |
| [`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md) | `tasks`, `on` | tasks, detached processes |

``` r

options(dMod2.cores = 2)
fits <- mstrust(obj, mu[free], fixed = gauge, sd = sigma[free], fits = 20,
                cores = c(fits = 4, conditions = 2))
profiles <- profile(obj, best, names(best), fixed = gauge, cores = 8)
```

## Tasks and job code

`distributedComputing(expr, tasks, on)` evaluates `expr` once per task
in a fresh R process; the value is the task result. `tasks` is a count,
a list or vector with one element per task, or a data frame with one row
per task. Inside `expr`:

- `task`: the element or row (as a list, `task$n`), `NULL` when `tasks`
  is a count,
- `taskId`: the position in `tasks`,
- `cores`: cores available to the task, e.g. for
  `mstrust(cores = cores)`, `parallel::mclapply(mc.cores = cores)` or
  `options(dMod2.cores = cores)`,
- `jobname`: the job name.

`input` names the session objects the code reads, functions included, by
default all objects of the calling environment. The process attaches
dMod2 and the R defaults; other packages are loaded by the code, and
session options such as `dMod2.cores` do not reach it.

Example: coverage of the 95% t interval for the mean of exponential
samples, one task per sample size, on the local machine:

``` r

interval <- function(y) t.test(y)$conf.int
sizes <- data.frame(n = c(5, 10, 20, 50))
job <- distributedComputing(
  {
    message("task ", taskId, ": samples of size ", task$n)
    covered <- replicate(2000, {
      ci <- interval(rexp(task$n))
      ci[1] < 1 && 1 < ci[2]
    })
    data.frame(n = task$n, coverage = mean(covered))
  },
  tasks = sizes, on = localNode(cores = 1), jobname = "coverage", input = "interval", seed = 1
)
#> coverage: 4 tasks started on localhost (1 core)
job
#> Job coverage: 4 tasks on localhost (1 core)
#> check(), wait(), get(combine = TRUE, all = FALSE), log(i), purge(local = FALSE), terminate()
```

The call returns once tasks have started. `job` is a list of functions:
`check()` reports task states, `wait()` polls until none runs,
[`get()`](https://rdrr.io/r/base/get.html) joins the results, here into
one data frame.

``` r

job$check()
#> coverage: 0 of 4 tasks done, 4 running
job$wait(interval = 1)
#> coverage: 0 of 4 tasks done, 4 running
#> coverage: 1 of 4 tasks done, 3 running
#> coverage: 2 of 4 tasks done, 2 running
#> coverage: 3 of 4 tasks done, 1 running
#> coverage: 4 of 4 tasks done
job$get()
#>    n coverage
#> 1  5   0.8760
#> 2 10   0.8960
#> 3 20   0.9285
#> 4 50   0.9445
```

Task `i` uses seed `seed + i`, independent of worker and order; without
`seed` one is drawn. A task that errors returns the `try-error`, the
others go on.

## Job folder layout

At submission the job files go into a temporary stage, copied to every
node as one `tar` stream over ssh into `<jobname>_folder/<k>` of the
home directory, `k` the node number. On `"localhost"` the folder is in
the working directory and no ssh is used. The local `<jobname>_folder`
also holds the manifest `job.rds`. After the local job:

``` r

list.files("coverage_folder", recursive = TRUE)
#>  [1] "1/job.R"           "1/slot_1.pid"      "1/slot_1.sh"      
#>  [4] "1/task_1.out"      "1/task_1.rds"      "1/task_2.out"     
#>  [7] "1/task_2.rds"      "1/task_3.out"      "1/task_3.rds"     
#> [10] "1/task_4.out"      "1/task_4.rds"      "1/tasks.rds"      
#> [13] "1/workspace.RData" "job.rds"
```

| file | content |
|----|----|
| `workspace.RData` | objects of `input`; `zstd`-compressed when available locally |
| `tasks.rds` | tasks as a list, one entry per data frame row |
| `job.R` | R script of a task |
| model files | sources, objects or shared objects, see [Model files and builds](#model-files-and-builds) |
| `build.sh`, `files.txt`, `build.log` | model build on the node, build output |
| `slot_<k>.sh`, `.pid`, `.run` | runner script of node `k`, process ID of the script and of the running task |
| `job.sh`, `env.sh` | SLURM only: job script from the template, job environment |
| `task_<i>.out`, `task_<i>.rds` | output and result of task `i` |

`job.R` is the same for every task, called as
`Rscript job.R <taskId> <cores>` with output in `task_<taskId>.out`:

``` r

writeLines(readLines("coverage_folder/1/job.R"))
#> library(dMod2)
#> load('workspace.RData')
#> for (so in list.files(pattern = paste0('[.]', substring(.Platform$dynlib.ext, 2L), '$'))) dyn.load(so)
#> taskId <- as.integer(commandArgs(TRUE)[1L])
#> cores <- as.integer(commandArgs(TRUE)[2L])
#> task <- readRDS('tasks.rds')[[taskId]]
#> jobname <- 'coverage'
#> set.seed(1 + taskId)
#> result <- try({
#> {
#>     message("task ", taskId, ": samples of size ", task$n)
#>     covered <- replicate(2000, {
#>         ci <- interval(rexp(task$n))
#>         ci[1] < 1 && 1 < ci[2]
#>     })
#>     data.frame(n = task$n, coverage = mean(covered))
#> }
#> })
#> saveRDS(result, paste0('task_', taskId, '.part'))
#> invisible(file.rename(paste0('task_', taskId, '.part'), paste0('task_', taskId, '.rds')))
```

It attaches dMod2, loads workspace and shared objects, reads the task,
sets the seed and evaluates the code in
[`try()`](https://rdrr.io/r/base/try.html), writing `task_<taskId>.rds`
for [`get()`](https://rdrr.io/r/base/get.html). Tasks share nothing but
files. The code runs in the job folder of the node; files written there
come back with `get(all = TRUE)`.

## Nodes and workers

`on` says where tasks run:

``` r

on <- localNode()                                        # the local machine, with the free cores
on <- localNode(cores = 8)
on <- sshNodes(c("server1", "server2"))                  # each node with the free cores there
on <- sshNodes(c("server1", "server2"), cores = 10)
on <- sshNodes(cores = c(server1 = 24, server2 = 20))
on <- sshNodes(c("localhost", "server1"))                # the local machine and one more
on <- slurmCluster("mycluster", cores = 8, walltime = "12:00:00")
```

A node is an ssh destination (`"host"`, `"user@host"`) with key
authentication, since no terminal takes a password. `"localhost"` runs
without ssh, in the session shell and working directory, with the R and
library paths of the session, so no second dMod2 install is needed;
[`localNode()`](https://dmodverse.github.io/dMod2/reference/localNode.md)
is `sshNodes("localhost")`.

Each node is one worker. Task `i` goes to worker `(i - 1) %% W + 1` of
`W`; a worker runs the assigned tasks in sequence, each with the `cores`
of the node. `tasks = NULL`, the default, gives one task per worker.

`cores` per task is one number, a vector named by node, or `NULL` for
the free cores at submission.
[`detectFreeCores()`](https://dmodverse.github.io/dMod2/reference/detectFreeCores.md)
counts cores minus the one-minute load average, on a local Windows
machine minus the processor load since the previous call, at least 1,
locally or over ssh; a node gets at most the free cores, with a warning
if given more.

``` r

detectFreeCores(c("localhost", "server1", "server2"))
```

On each node `slot_<k>.sh` calls `Rscript job.R` for each task of the
node in turn, started as `nohup setsid sh slot_<k>.sh &`, so it survives
the end of the ssh connection or the session. Tasks run with
`OMP_NUM_THREADS=1` and `MKL_NUM_THREADS=1`, so BLAS and OpenMP start
only threads the code asks for, and with `libs` first in the library
paths.

## Model files and builds

The model files of a job are those of the model objects in `input`
(`obsfn`, `parfn`, `prdfn`, objectives), found through their sources;
without any, all model files of `getOption("dMod2.outdir")`, else of the
working directory. `build` sets how they reach the nodes:

- `"compile"`: sources are copied and `build.sh` compiles them into one
  shared object against the cppDE installed on the node, with cppDE
  headers, SUNDIALS, KLU, the model `-D` switches and `-O2`.
  `buildCores` compilers run at once, by default the node cores, at
  most 16. Past the command-line limit objects go through a static
  archive, `buildBundle` sources per translation unit.
- `"link"`: object files are copied and linked, which needs a compatible
  compiler on the node.
- `"copy"`: shared objects are copied, which needs a compatible node.
- `"auto"`, the default: `"copy"` when all nodes are `"localhost"`, else
  `"compile"`.

The build runs once per node before the first task; `job.R` points every
model object of the workspace at the one shared object
([`modelname()`](https://dmodverse.github.io/dMod2/reference/modelname.md)).
On SLURM it runs once on the login node for all array tasks. Output goes
to `build.log`; a failed build stops the submission and shows the end of
the log. `libs` prepends library paths on the nodes for build and tasks,
expanded remotely, so `~` is the remote home.

The coverage job has no model files: no model object in `input`, no
shared object in the working directory.

## Job monitoring and recovery

| function | does |
|----|----|
| `check()` | counts done, running, failed; `TRUE` once none runs, warns on failures |
| `wait(interval = 10)` | polls until none runs, reports changes; `TRUE` if all done |
| `get(combine = TRUE, all = FALSE)` | returns the results |
| `log(i, lines = 20)` | last output lines of task `i` |
| `purge(local = FALSE)` | deletes the job folders on the nodes, on `"localhost"` the local `<jobname>_folder/<k>`, with `local = TRUE` also local folder and registry entry |
| `terminate()` | stops the tasks |

A task is done when the result file exists, running while the worker
runs (also while waiting), failed when the worker ended without a
result. [`get()`](https://rdrr.io/r/base/get.html) copies remote results
over ssh into the local `<jobname>_folder/<k>`, with `all = TRUE` all
files. It joins parlists with [`c()`](https://rdrr.io/r/base/c.html) and
data frames, parameter frames included, with
[`rbind()`](https://rdrr.io/r/base/cbind.html), which merges the two
sides of a profile; other results come as a list. Missing results are
dropped with a warning. `get(combine = FALSE)` returns one entry per
task, named like `tasks`: `NULL` if not done, the `try-error` if failed.
`terminate()` ends the process group of every slot, on SLURM the job by
`scancel`.

Output of the second task, then deletion:

``` r

job$log(2)
#> task 2: samples of size 10
job$purge(local = TRUE)
```

`distributedComputing(wait = TRUE)` waits, returns
[`get()`](https://rdrr.io/r/base/get.html) and deletes the job; a job
with failed tasks is kept for
[`log()`](https://rdrr.io/r/base/Log.html).

A job outlives the session. The manifest is in the local
`<jobname>_folder`, in the `<jobname>_folder` of every remote node until
`purge()`, and in the registry `tools::R_user_dir("dMod2", "data")`
until `purge(local = TRUE)`;
[`dMod2Jobs()`](https://dmodverse.github.io/dMod2/reference/dMod2Jobs.md)
lists registered jobs. `recover = TRUE` with a job name returns the job
functions without resubmitting, from the manifest in the working
directory, else the registry, else, with `on`, the remote node folders,
which recovers a job from another machine. A new job of the same name is
an error while the old one runs where the new one would replace the old
folders: in the same working directory or on the same remote nodes, also
from another machine.

``` r

dMod2Jobs()
job <- distributedComputing(jobname = "fits", recover = TRUE)
job <- distributedComputing(jobname = "fits", on = sshNodes(c("server1", "server2")),
                            recover = TRUE)
```

## SLURM clusters

[`slurmCluster()`](https://dmodverse.github.io/dMod2/reference/localNode.md)
names the login node, reached by ssh, and the resources per array task:
`cores`, `walltime`, `memPerCore` in GB. Workers are the `array` array
tasks, by default one per task; array task `j` runs tasks `j`,
`j + array`, `j + 2 * array`, … in sequence. The stage goes to
`~/<jobname>_folder/1` on the login node, where the template environment
is set up, the model files are built once, and `sbatch` submits from
that folder, so compute nodes must see the login-node home directory.
`partition`, `account`, `maxRunning` (array tasks at once) and `sbatch`
(further options such as `"--constraint=ib"`) override the `#SBATCH`
lines of the template. The start reports the job ID; `check()` reads
`squeue`, `terminate()` calls `scancel`.

[`slurmTemplate()`](https://dmodverse.github.io/dMod2/reference/slurmTemplate.md)
sets up the job script template once per login node, stored in
`tools::R_user_dir("dMod2", "config")`. A template is a cluster job
script with placeholders `{{ jobname }}`, `{{ ntasks }}` (array tasks),
`{{ cores }}`, `{{ memPerCore }}`, `{{ walltime }}`, and one line
`{{ run }}` that runs the tasks. Non-comment lines before `{{ run }}`
set up the environment, e.g. `module load`, also for the build on the
login node. A placeholder without value is an error, so `{{ walltime }}`
needs `walltime`. The shipped template for the bwForCluster Helix:

``` r

writeLines(readLines(system.file("slurm", "helix.sh", package = "dMod2")))
#> #!/bin/bash
#> # Job script template of distributedComputing() on slurmCluster("helix", ...) for
#> # the bwForCluster Helix, set up with slurmTemplate("helix"). The {{ ... }} fields
#> # are filled in, {{ run }} runs the tasks of one array task.
#> #SBATCH --job-name={{ jobname }}
#> #SBATCH --array=1-{{ ntasks }}
#> #SBATCH --partition=cpu-single
#> #SBATCH --nodes=1
#> #SBATCH --ntasks-per-node={{ cores }}
#> #SBATCH --mem-per-cpu={{ memPerCore }}gb
#> #SBATCH --time={{ walltime }}
#> #SBATCH --output=slurm_%A_%a.out
#> 
#> module load compiler/gnu/13.3
#> module load math/R
#> 
#> {{ run }}
```

``` r

slurmTemplate("mycluster", from = "~/mycluster.sh")   # once, from a job script, for editing
slurmTemplate("helix")                                 # the template dMod2 ships for Helix
slurmTemplate()                                        # the templates set up
cluster <- slurmCluster("helix", cores = 1, walltime = "12:00:00")
```

`slurmCluster(jobscript = )` takes another template for one job.

## Fits and profiles of the guide project

A cluster starts many small tasks sooner than one large one, so these
jobs use one core per task. The 200 fits of `fit.R`
([Optimization](https://dmodverse.github.io/dMod2/articles/guide-optimization.md))
as 200 one-fit tasks, each start drawn from the task seed:

``` r

free <- setdiff(names(mu), names(gauge))
fitsJob <- distributedComputing(
  mstrust(obj, mu[free], fixed = gauge, sd = sigma[free], fits = 1, iterlim = 1000, tolControl = tolerances),
  tasks = 200, on = cluster, jobname = "oppelt_fits", seed = 1,
  input = c("obj", "mu", "sigma", "free", "gauge", "tolerances"))
fitsJob$check()
fitsJob$wait(interval = 60)
fits <- as.parframe(fitsJob$get())
fitsJob$purge(local = TRUE)
```

The sources of `obj` compile once on the login node.
[`get()`](https://rdrr.io/r/base/get.html) joins the parlists,
[`as.parframe()`](https://dmodverse.github.io/dMod2/reference/as.parframe.md)
sorts by value. A later session, also on another day, takes up the job
by name:

``` r

fitsJob <- distributedComputing(jobname = "oppelt_fits", recover = TRUE)
```

The profiles of [Profile
likelihood](https://dmodverse.github.io/dMod2/articles/guide-profiles.md),
one side of one parameter per task:

``` r

sides <- expand.grid(par = names(best), side = c("left", "right"), stringsAsFactors = FALSE)
profilesJob <- distributedComputing(
  profile(obj, best, task$par, side = task$side, fixed = gauge,
          limits = c(lower = -7, upper = 7), stepControl = list(stop = "data"),
          algoControl = list(reoptimize = TRUE),
          optControl = list(iterlim = 100, tolControl = tolerances)),
  tasks = sides, on = cluster, jobname = "oppelt_profiles",
  input = c("obj", "best", "gauge", "tolerances"))
profilesJob$wait(interval = 60)
profiles <- profilesJob$get()
```

[`get()`](https://rdrr.io/r/base/get.html) binds both sides of each
parameter into one profile.

An [`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md)
multistart splits by `which`: every task draws all starts from the same
seed and fits those it names. The 200 starts of [Sparse
deviations](https://dmodverse.github.io/dMod2/articles/guide-sparse-deviations.md),
two per task:

``` r

compoundsJob <- distributedComputing({
  set.seed(3)
  msem(objCompounds + priorCompounds, random, init, fixed = gaugeReduced, fits = 200,
       sd = sdStart, keepFits = TRUE, which = task)
}, tasks = split(1:200, rep(1:100, each = 2)), on = cluster, jobname = "oppelt_compounds",
input = c("objCompounds", "priorCompounds", "random", "init", "gaugeReduced", "sdStart"))
fitCompounds <- compoundsJob$get()
```

[`get()`](https://rdrr.io/r/base/get.html) returns one
[`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md) fit, the
best, with table and fits of all starts. `helix.R` of the project runs
the other expensive guide scripts the same way.

Without SLURM, as many tasks run at once as there are nodes. One task
per node, the default, with ten fits per node core:

``` r

fitsJob <- distributedComputing(
  mstrust(obj, mu[free], fixed = gauge, sd = sigma[free], fits = 10 * cores, iterlim = 1000,
          tolControl = tolerances, cores = cores),
  on = sshNodes(c("server1", "server2")), jobname = "oppelt_fits", seed = 1,
  input = c("obj", "mu", "sigma", "free", "gauge", "tolerances"))
```

[Bile acid
transport](https://dmodverse.github.io/dMod2/articles/bileAcidTransport.md)
runs the multistarts as a local job, [Kok
2020](https://dmodverse.github.io/dMod2/articles/kokMSB2020.md) the 200
fits as a SLURM array job.

## System requirements

- A POSIX shell locally (on Windows the `sh` of Rtools) and on every
  remote node, bash on remote nodes for `build = "compile"` (default)
  and `"link"`.
- Remote nodes: ssh with key authentication, dMod2 and cppDE for the R
  on the node `PATH` or in `libs`, a compiler for `build = "compile"`.
- SLURM: an account on the login node, a template from
  [`slurmTemplate()`](https://dmodverse.github.io/dMod2/reference/slurmTemplate.md),
  a home directory visible to the compute nodes.
- `zstd` on every remote node if the local machine has it, since it then
  compresses the workspace.

Reference:
[`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md),
[`localNode()`](https://dmodverse.github.io/dMod2/reference/localNode.md),
[`sshNodes()`](https://dmodverse.github.io/dMod2/reference/localNode.md),
[`slurmCluster()`](https://dmodverse.github.io/dMod2/reference/localNode.md),
[`slurmTemplate()`](https://dmodverse.github.io/dMod2/reference/slurmTemplate.md),
[`detectFreeCores()`](https://dmodverse.github.io/dMod2/reference/detectFreeCores.md),
[`dMod2Jobs()`](https://dmodverse.github.io/dMod2/reference/dMod2Jobs.md),
[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md),
[`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md),
[`profile()`](https://dmodverse.github.io/dMod2/reference/profile.md)

[← 10. PEtab import and
export](https://dmodverse.github.io/dMod2/articles/guide-petab.md)
