# 11. Remote computing

Conditions of one evaluation run in threads, independent fits and
profiles in processes.
[`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md)
runs R code as tasks in detached processes, locally, over ssh or as a
SLURM array job. Only the local job below runs at page build.

## Threads and processes

`options(dMod2.cores = n)`, or `obj(pars, cores = n)` for one call,
solves the conditions of an evaluation in `n` OpenMP threads inside
cppDE; a single condition uses one.

[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md),
[`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md) and
[`profile()`](https://dmodverse.github.io/dMod2/reference/profile.md)
run fits as processes through `cores`. With one number, fits are forks
of the session with one thread each.
`cores = c(fits = 4, conditions = 2)` runs four fits with two threads
each in fresh R processes (a PSOCK cluster, always on Windows). Memory
grows with the process count; keep the product of both below the core
count.

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

`input` names the session objects the code reads, by default all objects
of the calling environment. The process attaches dMod2 and the R
defaults; session options do not reach it.

Example: coverage of the 95% t interval for exponential samples, one
task per sample size:

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

Task `i` uses seed `seed + i`. A task that errors returns the
`try-error`, the others go on.

## Job folder layout

The job files go to `<jobname>_folder/<k>` in the home directory of node
`k`, on `"localhost"` in the working directory. The local
`<jobname>_folder` also holds the manifest `job.rds`:

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

Every task runs `Rscript job.R <taskId> <cores>`:

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

The code runs in the job folder of the node; files written there come
back with `get(all = TRUE)`.

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
authentication. `"localhost"` runs without ssh, with the R and library
paths of the session;
[`localNode()`](https://dmodverse.github.io/dMod2/reference/localNode.md)
is `sshNodes("localhost")`.

Each node is one worker. Task `i` goes to worker `(i - 1) %% W + 1` of
`W`; a worker runs the assigned tasks in sequence with the `cores` of
the node. `tasks = NULL`, the default, gives one task per worker.
`cores` is one number, a vector named by node, or `NULL` for the free
cores, which
[`detectFreeCores()`](https://dmodverse.github.io/dMod2/reference/detectFreeCores.md)
counts as cores minus the load average.

``` r

detectFreeCores(c("localhost", "server1", "server2"))
```

Tasks survive the end of the ssh connection and the session, and run
with `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1` and `libs` first in the
library paths.

## Model files and builds

The model files of a job are those of the model objects in `input`;
without any, all model files of `getOption("dMod2.outdir")`, else of the
working directory. `build`:

- `"compile"`: sources are compiled on the node into one shared object
  against the cppDE installed there, with `buildCores` compilers
  (default the node cores, at most 16) and `buildBundle` sources per
  translation unit.
- `"link"`: object files are copied and linked, which needs a compatible
  compiler on the node.
- `"copy"`: shared objects are copied, which needs a compatible node.
- `"auto"`, the default: `"copy"` when all nodes are `"localhost"`, else
  `"compile"`.

The build runs once per node, on SLURM once on the login node. A failed
build stops the submission and shows the end of `build.log`. `libs`
prepends library paths on the nodes; `~` is the remote home.

## Job monitoring and recovery

| function | does |
|----|----|
| `check()` | counts done, running, failed; `TRUE` once none runs, warns on failures |
| `wait(interval = 10)` | polls until none runs, reports changes; `TRUE` if all done |
| `get(combine = TRUE, all = FALSE)` | returns the results |
| `log(i, lines = 20)` | last output lines of task `i` |
| `purge(local = FALSE)` | deletes the job folders on the nodes, on `"localhost"` the local `<jobname>_folder/<k>`, with `local = TRUE` also local folder and registry entry |
| `terminate()` | stops the tasks |

A task is failed when the worker ended without a result.
[`get()`](https://rdrr.io/r/base/get.html) copies remote results into
the local `<jobname>_folder/<k>`, with `all = TRUE` all files. It joins
parlists with [`c()`](https://rdrr.io/r/base/c.html) and data frames
with [`rbind()`](https://rdrr.io/r/base/cbind.html), which merges the
two sides of a profile; other results come as a list.
`get(combine = FALSE)` returns one entry per task: `NULL` if not done,
the `try-error` if failed.

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

A job outlives the session. The manifest is in the local and remote
`<jobname>_folder` and in the registry
`tools::R_user_dir("dMod2", "data")`;
[`dMod2Jobs()`](https://dmodverse.github.io/dMod2/reference/dMod2Jobs.md)
lists registered jobs. `recover = TRUE` with a job name returns the job
functions without resubmitting, from the working directory, the registry
or, with `on`, the remote nodes, also for a job from another machine. A
new job of the same name is an error while the old one runs in the same
folders.

``` r

dMod2Jobs()
job <- distributedComputing(jobname = "fits", recover = TRUE)
job <- distributedComputing(jobname = "fits", on = sshNodes(c("server1", "server2")),
                            recover = TRUE)
```

## SLURM clusters

[`slurmCluster()`](https://dmodverse.github.io/dMod2/reference/localNode.md)
names the login node and the resources per array task: `cores`,
`walltime`, `memPerCore` in GB. Workers are the `array` array tasks, by
default one per task. The job folder is `~/<jobname>_folder/1` on the
login node, where the model files are built once and `sbatch` submits,
so compute nodes must see that home directory. `partition`, `account`,
`maxRunning` and `sbatch` (further options) override the `#SBATCH` lines
of the template.

[`slurmTemplate()`](https://dmodverse.github.io/dMod2/reference/slurmTemplate.md)
sets up the job script template once per login node. A template is a job
script with placeholders `{{ jobname }}`, `{{ ntasks }}`, `{{ cores }}`,
`{{ memPerCore }}`, `{{ walltime }}` and one line `{{ run }}`; lines
before `{{ run }}` set up the environment, also for the build. The
shipped template for the bwForCluster Helix:

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

The 200 fits of
[Optimization](https://dmodverse.github.io/dMod2/articles/guide-optimization.md)
as 200 one-fit tasks:

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

A later session takes up the job by name:

``` r

fitsJob <- distributedComputing(jobname = "oppelt_fits", recover = TRUE)
```

The profiles of [Profile
likelihood](https://dmodverse.github.io/dMod2/articles/guide-profiles.md),
one side of one parameter per task;
[`get()`](https://rdrr.io/r/base/get.html) binds both sides into one
profile:

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

An [`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md)
multistart splits by `which`: every task draws all starts from the same
seed and fits those it names. The 200 starts of [Sparse
deviations](https://dmodverse.github.io/dMod2/articles/guide-sparse-deviations.md):

``` r

dcfJob <- distributedComputing({
  set.seed(3)
  msem(objSS, randomSS, c(initDCF, prob_compound = 0.2, tau_compound = 0.5),
       fixed = gaugeReduced, fits = 200, sd = 0.5, keepFits = TRUE, which = task)
}, tasks = as.list(1:200), on = cluster, jobname = "oppelt_dcf",
input = c("objSS", "randomSS", "initDCF", "gaugeReduced"))
fitDCF <- dcfJob$get()
```

[`get()`](https://rdrr.io/r/base/get.html) returns the best
[`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md) fit with
all starts. Over ssh, one task per node:

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
