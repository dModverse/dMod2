# Where distributedComputing() Runs Tasks

`localNode()` is the local machine, `sshNodes()` are machines reached by
`ssh` and `slurmCluster()` is a SLURM cluster reached through a login
node.
[`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md)
takes one of them as `on`.

## Usage

``` r
localNode(cores = NULL)

sshNodes(nodes = names(cores), cores = NULL)

slurmCluster(
  login,
  cores,
  array = NULL,
  walltime = NULL,
  memPerCore = 2,
  partition = NULL,
  account = NULL,
  maxRunning = NULL,
  sbatch = NULL,
  jobscript = NULL
)
```

## Arguments

- cores:

  The cores of every task. For `sshNodes()` one number for all nodes or
  a vector named by node. `NULL` (default) uses the free cores of each
  node, see
  [`detectFreeCores()`](https://dmodverse.github.io/dMod2/reference/detectFreeCores.md);
  a node gets at most the free cores it has. For `slurmCluster()` the
  cores of each array task.

- nodes:

  Character, the nodes, each once: `ssh` destinations such as `"host"`
  or `"user@host"` with key authentication, or `"localhost"` for the
  local machine. Defaults to the names of `cores`.

- login:

  Character, the login node, an `ssh` destination.

- array:

  The number of array tasks, at most one per task. `NULL` (default) is
  one per task.

- walltime:

  The maximum runtime of an array task, `"HH:MM:SS"`, filled into the
  job script. Default `NULL`, none.

- memPerCore:

  Memory per core in GB, filled into the job script. Default 2.

- partition, account:

  Character, the partition and the account of the job. Default `NULL`,
  those of the job script.

- maxRunning:

  The most array tasks running at once. Default `NULL`, no limit.

- sbatch:

  Character vector of further options of `sbatch`, such as
  `"--constraint=ib"`.

- jobscript:

  Path to a job script template, in place of the one
  [`slurmTemplate()`](https://dmodverse.github.io/dMod2/reference/slurmTemplate.md)
  set up for `login`. Default `NULL`.

## Value

An object of class `dMod2Nodes` or `dMod2Slurm`.

## Details

The tasks of a job run on workers: a node is one worker, and a SLURM job
has `array` workers, the array tasks of the job. The tasks are dealt to
the workers in turn, and a worker runs the tasks dealt to it one after
the other, each with `cores` cores. Without SLURM the processes run
detached and survive the R session. With SLURM the array job is
submitted on the login node through the job script template
[`slurmTemplate()`](https://dmodverse.github.io/dMod2/reference/slurmTemplate.md)
set up for it; `partition`, `account`, `maxRunning` and `sbatch` are
options of `sbatch`, which take precedence over the `#SBATCH` lines of
the template.

## See also

[`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md),
[`slurmTemplate()`](https://dmodverse.github.io/dMod2/reference/slurmTemplate.md)

## Examples

``` r
localNode(cores = 4)
#> localhost (4 cores)
sshNodes(c("server1", "server2"), cores = 10)
#> server1 (10 cores), server2 (10 cores)
sshNodes(cores = c(server1 = 24, server2 = 20))
#> server1 (24 cores), server2 (20 cores)
slurmCluster("mycluster", cores = 8, array = 32, walltime = "12:00:00")
#> SLURM cluster mycluster: 32 array tasks of 8 cores for 12:00:00
```
