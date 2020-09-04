### Code repository

The code in this folder is organized in several R scripts to obtain the data, run the analysis and produce the figures and tables. The scripts can be ran at once or in steps using GNU `make`. The only two requirements to reproduce this analysis are [git](https://git-scm.com) and [docker](https://www.docker.com).

The code uses two folders in the root directory of this repo; `data/` where the data would be downloaded and `article/` where the output would be saved.

The R scripts are

1. To get the data: `data-get.R`
2. To run the analysis: `target-run.R`
3. To produce the figures: `fig-ecdf.R` and `fig-variations.R`
4. To produce the tables: `tab-datasets.R`, `tab-heads.R` and `tab-cutoffs.R`

In addition, the `function-def.R` contains functions definitions that are used throughout.

To get and run the code, follow

```bash
git clone https://github.com/MahShaaban/ReScience-submission
cd ReScience-submission
make -f code/Makefile
```

The software environment where this code was ran and tested is publicly available as a [docker](https://www.docker.com) image.
To download and run the image, follow

```bash
docker pull bcmslab/target:20190927
docker run -it bcmslab/target:20190927 bash
```
NOTE: `article/` is currently a submodule to house the [Overleaf](https://www.overleaf.com/) private project. This should not cause any problems for running the code and reproducing the full analysis outputs; figures and tables.

This code is under GPL-3, see `LICENSE.md`.
