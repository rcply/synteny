
# Synteny

Data and scripts associated with the paper: *Sponges, ctenophores and the statistical significance of syntenies*

The program `rhduet.py` will do the hierarchical analysis, on the data in the `bed` and
`all_vs_all` (reciprocal best hits file) directories.

It will run as is for a comparison of *Capsaspora*, *Hormiphora*, *Ephydatia* and *Rhopilema*. Otherwise
inspect the arguments. 

Species are listed through `--sp_ids` and processed from *right* to *left*, so left should be the outgroup.

The option `--last` will restrict output to the final analysis round.