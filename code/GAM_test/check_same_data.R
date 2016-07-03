rm(list=ls())

## Check data similarity

steves <- read.csv("../../data/original/interp_short_allsystem_original.csv")
repros <- read.csv("../../data/reproduction/Repro_data_for_gam_test.csv", row.names=1)
name_table <- read.csv("../../data/reproduction/repro_steve_name_table.csv")
names(repros) <- name_table$steve_names[match(names(repros), name_table$repro_names)]

## Check Day.numbering
is.sorted(steves$Day.Number)
is.sorted(repros$Day.Number)

unique(diff(steves$Day.Number)) ## odd that there is one or more step of 6.8 days
table(round(diff(steves$Day.Number),2)) ## only one
steves$Day.Number[diff(steves$Day.Number)>6]

## tidy Steve's
from.steve <- gather(steves, Species, Abundance, 2:13)
names(from.steve) <- c("Day.number", "variable", "value")
from.steve$value_steve <- from.steve$value 
from.steve <- select(from.steve, -value)
## correcting an offset in the two datasets Day numbers
from.steve$Day.number <- from.steve$Day.number + 1.65 + 3.35

## tidy the repro data
from.repro <- gather(repros, Species, Abundance, 2:13)
names(from.repro) <- c("Day.number", "variable", "value")
from.repro$value_repro <- from.repro$value 
from.repro <- select(from.repro, -value)

## join the data, again first making days character, then back to numeric after
from.repro$Day.number <- as.character(from.repro$Day.number)
from.steve$Day.number <- as.character(from.steve$Day.number)
ff <- inner_join(from.repro, from.steve)
from.repro$Day.number <- as.numeric(from.repro$Day.number)
from.steve$Day.number <- as.numeric(from.steve$Day.number)
ff$Day.number <- as.numeric(ff$Day.number)


library(ggplot2)
ff$diff <- abs(ff$value_repro - ff$value_steve)
ff$diff_logs <- abs(log10(ff$value_repro) - log10(ff$value_steve))
ff$diff_colour <- ifelse(ff$diff_logs<(0.1), "less", "more")
ff$diff_colour[ff$diff==0] <- "less"
table(ff$diff_colour)

qplot(ff$diff)
range(ff$diff)

ggplot(filter(ff, variable=="TotalN"),
       aes(x=log10(value_steve), y=log10(value_repro), colour=diff_colour)) +
  geom_point() +
  facet_wrap(~variable, scales = "free") 

g1 <- ggplot(ff,
       aes(x=log10(value_steve), y=log10(value_repro), colour=diff_colour)) +
  geom_point()
g1
g1 +facet_wrap(~variable, scales = "free") + geom_abline(intercept=0, slope=1)


## and plot time series
g1 <- ggplot(ff, aes(x=as.numeric(Day.number), y=sqrt(value_repro))) +
  facet_wrap(~variable, ncol=2, scales="free_y") +
  geom_line(size=0.5, col="black") 
g2 <- geom_line(aes(x=Day.number, y=sqrt(value_steve)), colour="red")
g3 <- geom_point(aes(x=Day.number, y=sqrt(value_steve), colour=diff_colour), size=0.5)    
g1 + g2 + g3





