
# Synteny

Data and scripts associated with the paper: *Sponges, ctenophores and the statistical significance of syntenies*

The program `rhduet.py` will do the hierarchical analysis, on the data in the
`bed` and `all_vs_all` (reciprocal best hits file) directories. It requires the
python packages `numpy`,`statsmodels` and `scipy`. I used `python=3.10`, but 
I imagine it works with most recent v3 pythons.

It will run as is for a comparison of *Capsaspora*, *Hormiphora*, *Ephydatia* and *Rhopilema*. Otherwise
inspect the arguments. 

Species are listed through `--sp_ids` and processed from *right* to *left*, so left should be the outgroup.
e.g. `--sp_ids s_ros h_cal e_mue r_esc`

The option `--last` will restrict output to the final analysis round.
