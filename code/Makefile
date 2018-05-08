JLEXEC = julia

SRC = $(wildcard figure*.jl)
FIG = $(SRC:.jl=.pdf)

ALL: $(FIG)

%.pdf: %.jl
	$(JLEXEC) $<
