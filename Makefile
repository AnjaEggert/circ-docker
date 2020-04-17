RSCRIPT = Rscript --verbose $< &> $<.log

.PHONY: all

all: sims_proc2.rdata

sims.rdata: sims.R
	$(RSCRIPT)

sims_proc1.rdata: sims_proc1.R sims.rdata
	$(RSCRIPT)

sims_proc2.rdata: sims_proc2.R sims_proc1.rdata
	$(RSCRIPT)
