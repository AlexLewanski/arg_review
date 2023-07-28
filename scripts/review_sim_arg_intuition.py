##################################
##################################
### SIMULATIONS FOR ARG REVIEW ###
##################################
##################################

#Authors: Alex Lewanski, Mike Grundler, Gideon Bradburd
#Script written by Alex Lewanski


#####################
### SCRIPT SET-UP ###
#####################

### LIBRARIES ###
#import tskit
import msprime
import numpy as np
import pandas as pd


### PATHS, ETC ###
output_path = '/Users/alexlewanski/Documents/michigan_state/research/arg_review/figures/sim_material/output/'


### CUSTOM FUNCTIONS ###

def create_tree_height_df(tree_seq):
    tree_height_df = pd.DataFrame(columns=['tree_index', 'tree_height'])
    for TREE in tree_seq.trees():
        tree_height_df.loc[tree_height_df.shape[0]] = [TREE.index, TREE.time(TREE.root)]
    
    return tree_height_df

def create_tree_span_df(tree_seq):
    tree_span_df = pd.DataFrame(columns=['tree_index', 'left', 'right'])
    
    for TREE in tree_seq.trees():
        tree_span_df.loc[TREE.index] = [TREE.index, TREE.interval[0], TREE.interval[1]]
    
    return tree_span_df


def create_node_membership_df(tree_seq, sample_list = 'all'):
    if sample_list == 'all':
        sample_list = list(tree_seq.samples())

    node_membership_df = pd.DataFrame(columns=['sample', 'node_index', 'tree_index', 'node_time'])
    for SAMPLE in  sample_list:
        for TREE in full_arg_sim_simplify.trees():
            for NODE in TREE.nodes():
                if TREE.is_descendant(SAMPLE, NODE):
                    node_membership_df.loc[node_membership_df.shape[0]] = [SAMPLE, NODE, TREE.index, TREE.time(NODE)]
    
    return node_membership_df

def create_composition_df(tree_seq):
    node_composition_df = pd.DataFrame(columns=['node_index', 'tree_index', 'node_time'])
    for TREE in tree_seq.trees():
        for NODE in TREE.nodes():
            node_composition_df.loc[node_composition_df.shape[0]] = [NODE, TREE.index, TREE.time(NODE)]
    
    return node_composition_df



#################################################################
### PART I: SINGLE ARG SIMULATION TO LOOK AT GENERAL FEATURES ###
#################################################################

### SIMULATING AND SIMPLIFYING ARG ###
full_arg_sim = msprime.sim_ancestry(
    12,
    population_size = 100,
    recombination_rate = 0.00005, 
    sequence_length = 10000,
    record_full_arg = True, 
    random_seed = 11938)

full_arg_sim_simplify = full_arg_sim.simplify()


### PROCESSING SIMULATION OUTPUT ###
tree_height_df = create_tree_height_df(tree_seq = full_arg_sim_simplify)
tree_span_df = create_tree_span_df(tree_seq = full_arg_sim_simplify)
node_membership_df = create_node_membership_df(full_arg_sim_simplify, [0, 6, 11])
node_composition_df = create_composition_df(full_arg_sim_simplify)

breakpoint_simptree = [i for i in full_arg_sim_simplify.breakpoints()]
node_info_full_arg_df = pd.DataFrame({'time':full_arg_sim.nodes_time,
                                      'flags':full_arg_sim.nodes_flags})


### OUTPUTTING INFORMATION ###
tree_span_df.to_csv(output_path + 'tree_span_df.csv.gz', index = False, compression = "gzip")
node_membership_df.to_csv(output_path + 'node_membership_df.csv.gz', index = False, compression = "gzip")
node_composition_df.to_csv(output_path + 'node_composition_df.csv.gz', index = False, compression = "gzip")
tree_height_df.to_csv(output_path + 'tree_height_df.csv.gz', index = False, compression = "gzip")

full_arg_sim_simplify.write_nexus(output_path + 'sim_trees_review.nexus')

node_info_full_arg_df.to_csv(output_path + 'node_info_full_arg_df.csv.gz', index = False, compression = "gzip")




#########################################################################
### PART II: EXAMINING HOW DIFFERENT VARIABLES INFLUENCE ARG FEATURES ###
#########################################################################

### POPULATION SIZE SIMULATIONS ###
pop_size_tree_height_list = []
pop_size_tree_span_list = []
for POP_SIZE in np.arange(50, 1050, 50):
    for ITER in range(30):
        pop_size_sim = msprime.sim_ancestry(
            10,
            population_size = POP_SIZE,
            recombination_rate = 0.00003, 
            sequence_length = 10000,
            record_full_arg = False)
        pop_size_sim_simplify = pop_size_sim.simplify()

        pop_size_tree_height_df = create_tree_height_df(pop_size_sim_simplify)
        pop_size_tree_height_df['pop_size'] = POP_SIZE
        pop_size_tree_height_df['sim_index'] = ITER
        pop_size_tree_height_list.append(pop_size_tree_height_df)

        pop_size_tree_span_df = create_tree_span_df(pop_size_sim_simplify)
        pop_size_tree_span_df['pop_size'] = POP_SIZE
        pop_size_tree_span_df['sim_index'] = ITER
        pop_size_tree_span_list.append(pop_size_tree_span_df)

        print('pop size = ', POP_SIZE, '; iter: ', ITER)

pop_size_tree_height_df_combined = pd.concat(pop_size_tree_height_list)
pop_size_tree_span_combined = pd.concat(pop_size_tree_span_list)

pop_size_tree_height_df_combined.to_csv(output_path + 'pop_size_tree_height_df_combined.csv.gz', index = False, compression = "gzip")
pop_size_tree_span_combined.to_csv(output_path + 'pop_size_tree_span_combined.csv.gz', index = False, compression = "gzip")

len(np.arange(0, 0.1, 0.005))


### MIGRATION SIMULATIONS ###
mig_tree_height_list = []
mig_tree_span_list = []
for MIG_RATE in np.arange(0, 0.0001, 0.000005):
    for ITER in range(30):
        mig_demog = msprime.Demography()
        mig_demog.add_population(name = "pop1", initial_size = 500)
        mig_demog.add_population(name = "pop2", initial_size = 500)
        mig_demog.add_population(name = "ancestral_pop", initial_size = 500)
        mig_demog.add_population_split(time = 5000, derived = ["pop1", "pop2"], ancestral = "ancestral_pop")
        mig_demog.set_migration_rate(source = "pop1", dest = "pop2", rate = MIG_RATE)

        mig_sim = msprime.sim_ancestry(
            samples = {"pop1": 10},
            recombination_rate = 0.00003,
            sequence_length = 10000,
            record_full_arg = False,
            demography = mig_demog)
        
        mig_sim_simplify = mig_sim.simplify()

        mig_tree_height_df = create_tree_height_df(mig_sim_simplify)
        mig_tree_height_df['mig_rate'] = MIG_RATE
        mig_tree_height_df['sim_index'] = ITER
        mig_tree_height_list.append(mig_tree_height_df)

        mig_tree_span_df = create_tree_span_df(mig_sim_simplify)
        mig_tree_span_df['mig_rate'] = MIG_RATE
        mig_tree_span_df['sim_index'] = ITER
        mig_tree_span_list.append(mig_tree_span_df)

        print('mig rate = ', MIG_RATE, '; iter: ', ITER)

mig_tree_height_df_combined = pd.concat(mig_tree_height_list)
mig_tree_span_combined = pd.concat(mig_tree_span_list)

mig_tree_height_df_combined.to_csv(output_path + 'mig_tree_height_df_combined.csv.gz', index = False, compression = "gzip")
mig_tree_span_combined.to_csv(output_path + 'mig_tree_span_combined.csv.gz', index = False, compression = "gzip")


### SAMPLE SIZE SIMULATIONS ###
#***these results are not currently included in the manuscript
samp_size_tree_height_list = []
samp_size_tree_span_list = []
for SAMP_SIZE in np.arange(1, 59, 3):
    for ITER in range(30):
        samp_size_sim = msprime.sim_ancestry(
            SAMP_SIZE,
            population_size = 500,
            recombination_rate = 0.00005, 
            sequence_length = 10000,
            record_full_arg = False)
        samp_size_sim_simplify = samp_size_sim.simplify()

        samp_size_tree_height_df = create_tree_height_df(samp_size_sim_simplify)
        samp_size_tree_height_df['samp_size'] = SAMP_SIZE
        samp_size_tree_height_df['sim_index'] = ITER
        samp_size_tree_height_list.append(samp_size_tree_height_df)

        samp_size_tree_span_df = create_tree_span_df(samp_size_sim_simplify)
        samp_size_tree_span_df['samp_size'] = SAMP_SIZE
        samp_size_tree_span_df['sim_index'] = ITER
        samp_size_tree_span_list.append(samp_size_tree_span_df)

        print('sample size = ', SAMP_SIZE, '; iter: ', ITER)

samp_size_tree_height_combined = pd.concat(samp_size_tree_height_list)
samp_size_tree_span_combined = pd.concat(samp_size_tree_span_list)

samp_size_tree_height_combined.to_csv(output_path + 'samp_size_tree_height_combined.csv.gz', index = False, compression = "gzip")
samp_size_tree_span_combined.to_csv(output_path + 'samp_size_tree_span_combined.csv.gz', index = False, compression = "gzip")



#################################
### CODE NOT CURRENTLY IN USE ###
#################################
'''
tree_height_df = pd.DataFrame(columns=['tree_index', 'tree_height'])
for TREE in full_arg_sim_simplify.trees():
    tree_height_df.loc[tree_height_df.shape[0]] = [TREE.index, TREE.time(TREE.root)]

tree_span_df = pd.DataFrame(columns=['tree_index', 'left', 'right'])
for TREE in full_arg_sim_simplify.trees():
    #print(TREE.span)
    tree_span_df.loc[TREE.index] = [TREE.index, TREE.interval[0], TREE.interval[1]]

#for NODE in full_arg_sim_simplify.nodes():
#sample_array = full_arg_sim_simplify.samples()
node_membership_df = pd.DataFrame(columns=['sample', 'node_index', 'tree_index', 'node_time'])
for SAMPLE in [0, 6, 11]:
    for TREE in full_arg_sim_simplify.trees():
        for NODE in TREE.nodes():
            if TREE.is_descendant(SAMPLE, NODE):
                node_membership_df.loc[node_membership_df.shape[0]] = [SAMPLE, NODE, TREE.index, TREE.time(NODE)]
    print(SAMPLE)

node_composition_df = pd.DataFrame(columns=['node_index', 'tree_index', 'node_time'])
for TREE in full_arg_sim_simplify.trees():
    for NODE in TREE.nodes():
        node_composition_df.loc[node_composition_df.shape[0]] = [NODE, TREE.index, TREE.time(NODE)]


### TREE HEIGHT ###
tree_height_df = pd.DataFrame(columns=['tree_index', 'tree_height'])
for TREE in full_arg_sim_simplify.trees():
    tree_height_df.loc[tree_height_df.shape[0]] = [TREE.index, TREE.time(TREE.root)]

    
recombo_tree_height_list = []
recombo_tree_span_list = []
for RECOMBO in np.arange(0.00005, 0.001005, 0.00005):
    for ITER in range(50):
        recombo_sim = msprime.sim_ancestry(
            10,
            population_size = 1000,
            recombination_rate = RECOMBO, 
            sequence_length = 10000,
            record_full_arg = False)
        recombo_sim_simplify = full_arg_sim.simplify()

        recombo_tree_height_df = create_tree_height_df(recombo_sim_simplify)
        recombo_tree_height_df['recombo'] = RECOMBO
        recombo_tree_height_df['sim_index'] = ITER
        recombo_tree_height_list.append(recombo_tree_height_df)

        recombo_tree_span_df = create_tree_span_df(recombo_sim_simplify)
        recombo_tree_span_df['recombo'] = RECOMBO
        recombo_tree_span_df['sim_index'] = ITER
        recombo_tree_span_list.append(recombo_tree_span_df)

        print('recombo = ', RECOMBO, '; iter: ', ITER)

recombo_tree_span_combined = pd.concat(recombo_tree_span_list)
'''
