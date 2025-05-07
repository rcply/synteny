#!/usr/bin/env python

import sys,os
import random
import argparse
import statsmodels.stats.multitest as ssmulti
from scipy.stats import binom, poisson, hypergeom
from collections import defaultdict,Counter

###########################################

def main(args):

   randomize = 0

   sp_ids = args.sp_ids

   chr_lookup = defaultdict(dict)

   for sp_id in sp_ids:
      bed_file = args.bed_dir + sp_id + '.bed'
      with open(bed_file) as fh:
         for line in fh:
            line = line.rstrip()
            parts = line.split()
            chr = parts[0]
            gene_id = parts[3]
            chr_lookup[sp_id][gene_id] = chr

   ingroup  = len(sp_ids) - 1
   outgroup = ingroup - 1

   valid_set = None

   doing_random = random_counter = 0
   while(outgroup >= 0 or doing_random):
      in_sp_ids = sp_ids[ingroup:]
      out_sp_id = sp_ids[outgroup]
      print(in_sp_ids,'vs',out_sp_id)
      all_species = [out_sp_id] + in_sp_ids
      orth_fn = write_orthologs(all_species,args.rbh_dir)
      orth_groups = read_orth_groups(all_species,orth_fn,chr_lookup)
      # permutations here
      if doing_random:
         print("doing randomization:",random_counter,out_sp_id,file=sys.stderr)
         chr_lookup = shuffle(out_sp_id,chr_lookup)
      equivalences,gene_locations = get_equivalences(orth_groups,chr_lookup)
      if args.show_loc and outgroup == 0:
         print_ortholog_chrs(gene_locations)
      counts = count_equivalences(equivalences)
      p_vals,out_set = get_p_vals(equivalences,counts,valid_set)

      print("n_tests: ",len(p_vals))
      p_adjs = ssmulti.multipletests(p_vals, alpha=float(args.alpha), method="fdr_bh")
#      p_adjs =  ssmulti.fdrcorrection_twostage(p_vals,alpha=0.05,method="bky")
#      p_adjs = ssmulti.multipletests(p_vals, alpha=0.1, method="fdr_bh")
      if outgroup > 0 and args.show_last == True:
         no_print = 1
      else:
         no_print = 0
      if outgroup > 0:
         valid_set = update_valid_set(equivalences,out_set,counts,p_vals,p_adjs,no_print,random_counter)
      else:
         update_valid_set(equivalences,out_set,counts,p_vals,p_adjs,no_print,random_counter)

      # deciding what to do next...

      if outgroup == 0 and doing_random == 0 and args.do_rand:
         doing_random = 1
         random_counter = 1
      elif doing_random:
         if random_counter == int(args.do_rand):
            sys.exit()
         random_counter += 1
      else:
         outgroup -= 1
         ingroup -= 1

###########################################

def shuffle(shuffle_id,chr_lookup):

   shuff_genes = list(chr_lookup[shuffle_id].keys())
   random.shuffle(shuff_genes)
   chr_lookup[shuffle_id] = dict(zip(shuff_genes,chr_lookup[shuffle_id].values()))
   return chr_lookup

###########################################

def update_valid_set(equivalences,out_set,counts,p_vals,p_adjs,no_print,random_counter):

   new_set = set()
   for i, o in enumerate(out_set):
       chr_list = [ gene[1] for gene in o ]
       chr_list = tuple(chr_list)
       if p_adjs[0][i] == True:
          new_set.add(chr_list)
       if no_print: continue
       total_ig = counts[chr_list[1:]]
       total_og = counts[chr_list[0]]
       total_me = equivalences[o]
       chr_str = ' '.join([ str(gene[1]).rjust(6) for gene in o ])
       sp_str = ' '.join([ gene[0] for gene in o ])
       if random_counter:
          tag = 'R' + str(random_counter) + ':' + str(i) + ':'
          print(tag,end=' ')
       print("{:8.2g}".format(p_adjs[1][i]),str(p_adjs[0][i]).rjust(6),
             "{:8.2g} {:3d} {:3d} {:3d}".format(p_vals[i], total_og, total_me, total_ig),chr_str,sp_str)
   return new_set

###########################################

def get_p_vals(equivalences,counts,valid_set):

   p_vals = []
   out_set = []

   n_ogs = equivalences.total()

   for tup, n_obs in equivalences.most_common():
      chr_list = [ gene[1] for gene in tup ]
      chr_list = tuple(chr_list)
      og_tag = chr_list[0]
      ig_tag = tuple(chr_list[1:])
      if not valid_set or (valid_set and (ig_tag in valid_set)):
         ph_val = hypergeom.sf(n_obs - 1, n_ogs, counts[og_tag], counts[ig_tag])
         p_vals.append(ph_val)
         out_set.append(tup)

   return p_vals,out_set

###########################################

def count_equivalences(equivalences):

    counts = Counter()

    for tup, n_obs in equivalences.most_common():
       chr_list = [ gene[1] for gene in tup ]
       chr_list = tuple(chr_list)
       counts[chr_list[0]] += n_obs
       counts[chr_list[1:]] += n_obs

    return counts

###########################################

def print_ortholog_chrs(locations):

   for orths in locations:
      chrs = [ x[1].rjust(6) for x in orths ]
      genes = [ x[0] for x in orths ]
      print(' '.join(chrs),' '.join(genes))
   return

###########################################

def get_equivalences(orth_groups,chr_lookup):

   observed_equivalences = Counter()
   all_locations = []

   for orth_group in orth_groups:
      tuple_list = []
      location_list = []
      for gene in orth_group:
         sp_id = gene[0]
         gene_id = gene[1]
         chr_id = chr_lookup[sp_id][gene_id]
         tuple_list.append((sp_id,chr_id))
         location_list.append((gene_id,chr_id))
      observed_equivalences[tuple(tuple_list)] += 1
      all_locations.append(location_list)

   return observed_equivalences,all_locations

###########################################

def read_orth_groups(sp_ids,orths_file,chr_lookup):

   with open(orths_file) as orths:

      orth_groups = []
      for line in orths:
         line = line.rstrip()
         parts = line.split()
         orth_counter = 0
         orth_group = []
         for orth_id in parts:
            sp_id,gene_id = orth_id.split('|')
            if gene_id not in chr_lookup[sp_id]:
               pass
            else:
               orth_group.append((sp_id,gene_id))
               orth_counter += 1
         if orth_counter == len(sp_ids):
            orth_group = sorted(orth_group,key=lambda x: sp_ids.index(x[0]))
            orth_groups.append(orth_group)

   return orth_groups

###########################################

def write_orthologs(sp_list,rbh_dir):

   out_fn = '_'.join(sp_list) + '_orthologs.txt'
   out_fh = open(out_fn,'w')

   rbh_by_gene = defaultdict(set)

   species = set()

   for sp_a in sp_list:
      for sp_b in sp_list:
         if sp_a == sp_b: continue
         rbh_file = rbh_dir + sp_a + '_nr_vs_' + sp_b + '_nr_rbh.txt'
         if not os.path.isfile(rbh_file):
            continue
         species.add(sp_a)
         species.add(sp_b)
         with open(rbh_file) as fh:
            for line in fh:
               line = line.rstrip()
               parts = line.split()
               gene_a = sp_a + '|' + parts[0]
               gene_b = sp_b + '|' + parts[1]
               rbh_by_gene[gene_a].add(gene_a)
               rbh_by_gene[gene_a].add(gene_b)
               rbh_by_gene[gene_b].add(gene_b)
               rbh_by_gene[gene_b].add(gene_a)

   n_species = len(species)

   done = dict()

   for gene in rbh_by_gene:
      gene_set = rbh_by_gene[gene]
      if len(gene_set) != n_species: continue
      fail = 0
      new_set = set()
      for member in gene_set:
         sp_set = { x.split('|')[0] for x in rbh_by_gene[member] }
         if len(sp_set) != n_species: fail = 1
         new_set = new_set.union(rbh_by_gene[member])
         if len(new_set) != n_species: fail = 1
      if fail: continue
      new_set = sorted(new_set)
      if(str(new_set) not in done):
         print(' '.join(list(new_set)),file=out_fh)
         done[str(new_set)] = 1

   return out_fn

###########################################

if __name__ == "__main__":

   parser = argparse.ArgumentParser(description = '')
   parser.add_argument('--bed_dir',action='store',dest='bed_dir',
                       help='dir location of bed files',default='../renamed_bed/')
   parser.add_argument('--rbh_dir',action='store',dest='rbh_dir',
                       help='dir location of rbh files',default='../all_vs_all/')
   parser.add_argument('--alpha',action='store',dest='alpha',default=0.05,
                       help='alpha')
   parser.add_argument('--rand',action='store',dest='do_rand',default=0,
                       help='do 100 randomizations')
   parser.add_argument('--last',action='store_true',dest='show_last',
                       help='only show final result')
   parser.add_argument('--sp_ids',nargs='+',dest='sp_ids',
                       help='list of species ids',
                       default=['c_owc','h_cal','e_mue','r_esc'])
   parser.add_argument('--show_locations',action='store_true',dest='show_loc',
                       help='show ortholog locations')

   args = parser.parse_args()

   main(args)

