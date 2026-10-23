if (.Platform$OS.type == "unix") {
  # Two tasks on the local machine, in a temporary working directory
  old <- setwd(tempdir())
  M <- matrix(rnorm(1e2), 10, 10)
  job <- distributedComputing(det(M * task), tasks = list(a = 1, b = 2),
                              input = "M", jobname = "example")
  while (!job$check()) Sys.sleep(1)
  print(job$get())
  job$purge(local = TRUE)
  setwd(old)
}

\dontrun{
# For an objective obj and a centre of the starts: multistarts on two machines
# reached by ssh, one task per machine with ten fits per core
jobFitting <- distributedComputing(
  mstrust(obj, center, fits = 10 * cores, cores = cores, sd = 4),
  machine = c(server1 = 8, server2 = 16), jobname = "fits", compile = TRUE
)
while (!jobFitting$check()) Sys.sleep(10)
parlists <- Filter(function(f) inherits(f, "parlist"), jobFitting$get())
fits     <- as.parframe(do.call(c, parlists))
bestfit  <- as.parvec(fits)
jobFitting$purge()

# Profiles of the best fit, four parameters per task, on the free cores of
# both machines
jobProfiling <- distributedComputing(
  profile(obj, bestfit, whichPar = task, limits = c(-5, 5), cores = cores),
  tasks = split(names(bestfit), ceiling(seq_along(bestfit) / 4)),
  machine = c("server1", "server2"), jobname = "profiles"
)
while (!jobProfiling$check()) Sys.sleep(10)
profiles <- do.call(rbind, jobProfiling$get())

# A SLURM array job of 10 tasks with 20 cores each, on the login node
# "mycluster". Its job script template is set up once, from a job script of
# the cluster with the placeholders of ?slurmTemplate.
slurmTemplate("mycluster", from = "~/mycluster.sh")
job <- distributedComputing(
  mstrust(obj, center, fits = 20, cores = cores, sd = 4),
  tasks = 10, machine = c(mycluster = 20), slurm = TRUE,
  walltime = "02:00:00", jobname = "slurmfits", compile = TRUE
)

# The functions of a job submitted earlier from this working directory
jobProfiling <- distributedComputing(NULL, jobname = "profiles", recover = TRUE)
}
