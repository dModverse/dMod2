# Job script templates of SLURM clusters

Sets up the job script template of a SLURM cluster, once per login node.
The templates are kept in `tools::R_user_dir("dMod2", "config")` and
stay across sessions;
[`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md)
on a
[`slurmCluster()`](https://dmodverse.github.io/dMod2/reference/localNode.md)
takes the one named like the login node of the cluster.

## Usage

``` r
slurmTemplate(
  login = NULL,
  from = login,
  edit = interactive(),
  overwrite = FALSE
)
```

## Arguments

- login:

  The login node, as in
  [`slurmCluster()`](https://dmodverse.github.io/dMod2/reference/localNode.md).
  `NULL` (default) lists the templates set up.

- from:

  The template to start from: the name of one dMod2 ships, such as
  `"helix"`, or the path of a job script. Defaults to the shipped
  template named `login`.

- edit:

  Logical, opens the template in the editor. Default
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

- overwrite:

  Logical. `TRUE` replaces a template set up before; with `FALSE`
  (default) it is kept.

## Value

The path of the template, invisibly. Without `login` the paths of the
templates set up, named by login node.

## Details

A template is a job script with the placeholders `{{ jobname }}`,
`{{ ntasks }}` (the number of array tasks), `{{ cores }}`,
`{{ memPerCore }}`, `{{ walltime }}` and one line `{{ run }}`, which
runs the tasks of the array task `SLURM_ARRAY_TASK_ID`. The lines before
`{{ run }}` that are not comments set up the environment, for example
with `module load`; the build on the login node runs in it as well. A
job script of the cluster becomes a template by these placeholders.
dMod2 ships one for the bwForCluster Helix, `"helix"`.

## See also

[`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md),
[`slurmCluster()`](https://dmodverse.github.io/dMod2/reference/localNode.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# once per login node, from a job script of the cluster
slurmTemplate("mycluster", from = "~/mycluster.sh")
# the template dMod2 ships for the bwForCluster Helix
slurmTemplate("helix")
slurmTemplate()
} # }
```
