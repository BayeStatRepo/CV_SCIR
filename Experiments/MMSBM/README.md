# MMSBM experiments

### Simulation

To replicate the simulation results shown in Figure 4a of the paper, run the script `R/Simulation_study_MMSBM_KL.R`.

### Email-Eu-Core Network Dataset

The popular Email-Eu-Core network dataset is publicly available at [SNAP](http://snap.stanford.edu/data/). The network was generated using email data from a large European research institution (Leskovec et al., 2007).  

- An edge `(i, j)` exists if person `i` sent at least one email to person `j`.  
- Emails represent only communication between institution members (the “core”), for a total of **N = 1,005 members** represented as nodes.
- The dataset does not include emails to or from individuals outside the institution.  
- Ground-truth community memberships are provided: each individual belongs to exactly one of **K = 42 departments** at the research institute.  

The dataset is provided in the `Data/` directory.

## Reference

J. Leskovec, J. Kleinberg, and C. Faloutsos. *Graph evolution: Densification and shrinking diameters.* ACM Transactions on Knowledge Discovery from Data, 1(1):2es, Mar. 2007. doi: [10.1145/1217299.1217301](https://doi.org/10.1145/1217299.1217301)
