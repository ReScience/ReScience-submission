
# Replication of Hinton & Nowlan (1987) "How Learning Can Guide Evolution"

This is a repository with the code and article submitted to the [Re**Science** journal](https://rescience.github.io). 

The reference for the original article (APA style):  
[Hinton, G. E., & Nowlan, S. J. (1987). How learning can guide evolution. Complex systems, 1(3), 495-502.](http://www.complex-systems.com/pdf/01-3-6.pdf)

The structure of the repository is the following. Folder `article` contains a brief description of my replication of the evolutionary simulations from the original article (see `article/stojic2017.pdf`).

Folder `code` contains main script `simulation.R`, written in R that runs the simulation, produces the summary data in CSV format that is saved in `data` folder. Script `utils.R` loads required packages and defines additional useful functions. Data in `data` folder is used by `genFigure2.R` in the `code` folder to produce the figures used in the replication article. Folder `notebook` does not contain anything.


## How to run the simulations yourself?

The code has been developed on Linux operating system. To execute the code you will need to install [R](https://www.r-project.org/) and couple of R packages. More details on operating system and R is available in `code/README.md`.

Note:  
1. Simulations take some time - on a 16 core computer it takes about 8h. For this reason, the data is already provided in `data` directory.   
2. If you execute `code/simulation.R` it will automatically try to install needed versions of the packages and potentially overwrite your version if you have some of them installed. Code might very well work with your version but I cannot guarantee it.  