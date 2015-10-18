---
Title: "Least-cost modelling on irregular landscape graphs"
Author:
  - name: Joseph Stachelek
    affiliation: 1
Address:
  - code:    1
    address: South Florida Water Management District, West Palm Beach, Florida, USA
Contact:
  - jstachel@sfwmd.gov
Editor:
  - Name Surname
Reviewer:
  - Name Surname
  - Name Surname
Publication:
  received:  Sep,  1, 2015
  accepted:  Sep, 1, 2015
  published: Sep, 1, 2015
  volume:    "**1**"
  issue:     "**1**"
  date:      Sep 2015
Repository:
  article:   "http://github.com/rescience/rescience-submission/article"
  code:      "http://github.com/rescience/rescience-submission/code"
  data:      
  notebook:  
Reproduction:
  - Least-cost modelling on irregular landscape graphs, T. Etherington, Landscape Ecology, 2012
Bibliography:
  article.bib

---

# Introduction

We propose a reference implementation of [@etherington2012least] that introduces a method for generating accumulated cost surfaces using irregular landscape graphs. According to the original article, irregular landscape graphs allow for faster processing speeds and avoid directional bias relative to regular landscape graphs. The original implementation was made in Python whose sources are available upon request to the author of the original article. The reference implementation we propose has been coded in R because of the strength of existing libraries for generating accumulated cost surfaces using regular landscape graphs [@gdistance]. 

# Methods

We used the description of the model as well as the source code of the original implementation (requested from the author) as the basis for the following reference implementation. The original article contrasted the processsing speeds and directional bias of accumulated cost surfaces generated with irregular landscape graphs relative to regular landscape graphs. We attempted to follow the structure, style, and order-of-operations of the original with a few exceptions. For example, we use the same underlying algorithm for computing the initial Delaunay triangulations [@barber1996quickhull] that form the basis of irregular landscape graph construction. One notable difference in the reference implementation relative to the original is that we have used matrix operations from the gdistance package [@gdistance] rather than nested loops to contruct regular landscape graphs. 
# Results

Results should be compared with original results and you have to explain why
you think they are the same or why they may differ (qualitative result vs
quantitative result). Note that it is not necessary to redo all the original
analysis of the results.


# Conclusion

Conclusion, at the very minimum, should indicate very clearly if you were able
to replicate original results. If it was not possible but you found the reason
why (error in the original results), you should exlain it.


Heading 1                          Heading 2
---------- ----------- ----------- ----------- ----------- -----------
cell1 row1 cell2 row 1 cell3 row 1 cell4 row 1 cell5 row 1 cell6 row 1
cell1 row2 cell2 row 2 cell3 row 2 cell4 row 2 cell5 row 2 cell6 row 2
cell1 row3 cell2 row 3 cell3 row 3 cell4 row 3 cell5 row 3 cell6 row 3
---------- ----------- ----------- ----------- ----------- -----------

Table: Table caption {#tbl:table}

A reference to table @tbl:table.
A reference to figure @fig:logo.
A reference to equation @eq:1.
A reference to citation @markdown.

![Figure caption](rescience-logo.pdf) {#fig:logo}

$$ A = \sqrt{\frac{B}{C}} $$ {#eq:1}


# References
