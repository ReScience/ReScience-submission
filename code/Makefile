#!/bin/bash


all: log \
	log/data-get.R.Rout \
	log/target-run.R.Rout \
	article/figures/ecdf.png \
	article/figures/variations.png \
	article/tables/heads.tex \
	article/tables/cutoffs.tex \
	article/tables/datasets.tex
log: 
	test ! -d log/ && mkdir log || exit 0
log/data-get.R.Rout: code/data-get.R
	R CMD BATCH --vanilla $< log/$(<F).Rout
log/target-run.R.Rout: code/target-run.R \
	code/function-def.R \
	data/beta_input/ESR1_diff_expr.xls \
	data/beta_input/Tet1_diff_expr.xls \
	data/beta_input/AR_diff_expr.xls \
	data/beta_input/349_peaks.bed \
	data/beta_input/5795_peaks.bed \
	data/beta_input/3656_peaks.bed \
	data/beta_input/mm9.refseq \
	data/beta_input/hg19.refseq
	R CMD BATCH --vanilla $< log/$(<F).Rout
article/figures/ecdf.png: code/fig-ecdf.R \
	code/function-def.R \
	data/target_output/variants_Tet1_5795.rds \
	data/target_output/variants_ESR1_349.rds \
	data/target_output/variants_AR_3656.rds
	R CMD BATCH --vanilla $< log/$(<F).Rout
article/figures/variations.png: code/fig-variations.R \
	code/function-def.R \
	data/target_output/basic_Tet1_5795.rds \
	data/target_output/basic_ESR1_349.rds \
	data/target_output/basic_AR_3656.rds
	R CMD BATCH --vanilla $< log/$(<F).Rout
article/tables/heads.tex: code/tab-heads.R \
	code/function-def.R \
	data/target_output/basic_Tet1_5795.rds \
	data/target_output/basic_ESR1_349.rds \
	data/target_output/basic_AR_3656.rds
	R CMD BATCH --vanilla $< log/$(<F).Rout
article/tables/cutoffs.tex: code/tab-cutoffs.R \
	code/function-def.R \
	data/target_output/basic_Tet1_5795.rds \
	data/target_output/basic_ESR1_349.rds \
	data/target_output/basic_AR_3656.rds
	R CMD BATCH --vanilla $< log/$(<F).Rout
article/tables/datasets.tex: code/tab-datasets.R 
	R CMD BATCH --vanilla $< log/$(<F).Rout
