import mdtraj as md
import numpy as np

t = md.load("traj_293.xtc",top="Hg.pdb")

pairs_Hg_Hg = t.topology.select_pairs('name Hg','name Hg')
rdf_list = md.compute_rdf(t, pairs = pairs_Hg_Hg, r_range = (0.0,2.00000/2), bin_width = 0.002, periodic=True)
for r, rdf in zip(rdf_list[0], rdf_list[1]):
    print(r*10, rdf)
