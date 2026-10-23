if (.Platform$OS.type == "unix") {
  # Two tasks on the local machine, in a temporary working directory
  old <- setwd(tempdir())
  M <- matrix(rnorm(1e2), 10, 10)
  job <- distributedComputing(det(M * task), tasks = list(a = 1, b = 2),
                              input = "M", jobname = "example")
  job$wait(interval = 1)
  print(job$get())
  job$purge(local = TRUE)
  setwd(old)
}

\dontrun{
# For an objective obj and a centre of the starts: multistarts on two nodes
# reached by ssh, one task per node with ten fits per core
nodes <- sshNodes(c("server1", "server2"), cores = 16)
jobFitting <- distributedComputing(
  mstrust(obj, center, fits = 10 * cores, cores = cores, sd = 4),
  on = nodes, jobname = "fits"
)
jobFitting$wait()
fits    <- as.parframe(jobFitting$get())
bestfit <- as.parvec(fits)
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
# through the job script template of its login node "mycluster", set up once
# from a job script of the cluster with the placeholders of ?slurmTemplate.
slurmTemplate("mycluster", from = "~/mycluster.sh")
cluster <- slurmCluster("mycluster", cores = 10, walltime = "02:00:00")
job <- distributedComputing(
  mstrust(obj, center, fits = cores, cores = cores, sd = 4),
  tasks = 100, on = cluster, jobname = "slurmfits"
)

# The jobs submitted and not purged, and the functions of one of them
dMod2Jobs()
job <- distributedComputing(NULL, jobname = "slurmfits", recover = TRUE)
}
