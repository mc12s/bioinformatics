#/bin/Rscript
#usage: Rscript ConvergenceStatistics T_40000_2500/1_2_200_loci/r0/512/u_.0015/
#This script calculates the autocorrellation, effective sample size, and gelman rubin dignostic for each MCMC run

library(coda)

args <- commandArgs(trailingOnly = TRUE)
print(args)

infolder <- args[1]

for(state in c("200K","1M","5M","20M")){

    gs_total = 0
    ess_total = 0

    for(rep in seq(1:10)){

        a <- as.vector(read.table(paste(infolder, "/rep", rep, "/beast_", state, "_msout_1/posterior.txt", sep="")))
        b <- as.vector(read.table(paste(infolder, "/rep", rep, "/beast_", state, "_msout_2/posterior.txt", sep="")))

        write(autocorr(as.mcmc(a), 1:100), file=paste(infolder, "/rep", rep, "/beast_", state, "_msout_1/autocorr.txt", sep=""), ncolumns=1)
        write(autocorr(as.mcmc(b), 1:100), file=paste(infolder, "/rep", rep, "/beast_", state, "_msout_2/autocorr.txt", sep=""), ncolumns=1)

        ess_total = ess_total + effectiveSize(as.mcmc(a))
        ess_total = ess_total + effectiveSize(as.mcmc(b))

        write(effectiveSize(as.mcmc(a)), file=paste(infolder, "/rep", rep, "/beast_", state, "_msout_1/ess.txt", sep=""))
        write(effectiveSize(as.mcmc(b)), file=paste(infolder, "/rep", rep, "/beast_", state, "_msout_2/ess.txt", sep=""))
        pmc <- mcmc.list(as.mcmc(a$V1),as.mcmc(b$V1))

        gs <- gelman.diag(pmc)
        gs_total = gs_total + gs$psrf[2]
        write(gs$psrf[2], file=paste(infolder, "/rep", rep, "/beast_", state, "_msout_1/gs.txt", sep=""))
        # greater than 1.1 is non-convergence

    }

    ess_avg = ess_total/20
    write(ess_avg, file=paste(infolder, "/ess_avg_",state,".txt", sep=""))
    gs_avg = gs_total/10
    write(gs_avg, file=paste(infolder, "/grs_avg",state,".txt", sep=""))
}

