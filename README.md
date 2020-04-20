# circ-docker

circ docker Jupyter

run following

    docker build . -t circ-im
    docker run --name="circ-dock" -p 8888:8888 -v .:/home/jovyan/analysis circ-im

copy URL into browser. Jupyter notebook opens.
go to new -> Terminal
change to folder <analysis>

    cd ~/analysis

Install R package from local tar.gz

    Rscript -e 'install.packages("HarmonicRegression_1.9999.tar.gz")'

You can now run the R scripts manually in the terminal, one after the other:

    Rscript sims.R &> sims.log
    Rscript sims_proc1.R &> sims_proc1.log
    Rscript sims_proc2.R &> sims_proc2.log

It works
the Makefile does not work properly

    make
