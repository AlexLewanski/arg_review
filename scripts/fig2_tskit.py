import tskit
import numpy as np
from IPython.display import SVG,display

#####################################################
### CREATING THE TREE SEQUENCE IN TSKIT FOR FIG 2 ###
#####################################################

### EMPTY TABLES WITH ARBITRARY LENGTH OF 100 ###
tables = tskit.TableCollection(sequence_length=100)

### NODE TABLE ###
node_table = tables.nodes

#samples nodes
node_table.add_row(flags=tskit.NODE_IS_SAMPLE)  # Node 0 (defaults to time 0)
node_table.add_row(flags=tskit.NODE_IS_SAMPLE)  # Node 1 (defaults to time 0)
node_table.add_row(flags=tskit.NODE_IS_SAMPLE)  # Node 2 (defaults to time 0)
node_table.add_row(flags=tskit.NODE_IS_SAMPLE)  # Node 3 (defaults to time 0)

#non-sample nodes
node_table.add_row(time=5)  # Node 4 (not a sample)
node_table.add_row(time=10)  # Node 5 (not a sample)
node_table.add_row(time=15)  # Node 6 (not a sample)
node_table.add_row(time=20)  # Node 7 (not a sample)
node_table.add_row(time=25)  # Node 8 (not a sample)


### EDGE TABLE ###
edge_table = tables.edges
edge_table.set_columns(
    left=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 30.0, 30.0, 30.0, 30.0, 60.0, 60.0]),
    right=np.array([60.0, 30.0, 30.0, 60.0, 30.0, 30.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]),
    parent=np.array([8, 4, 4, 8, 6, 6, 6, 6, 5, 5, 7, 7], dtype=np.int32),  # References IDs in the node table
    child=np.array([0, 1, 2, 6, 3, 4, 1, 5, 2, 3, 0, 6], dtype=np.int32),  # References IDs in the node table
)


### SITES TABLE ###
tables.sites.add_row(position=8, ancestral_state='G')
tables.sites.add_row(position=16, ancestral_state='T')
tables.sites.add_row(position=24, ancestral_state='G')
tables.sites.add_row(position=32, ancestral_state='G')
tables.sites.add_row(position=40, ancestral_state='T')
tables.sites.add_row(position=48, ancestral_state='C')
tables.sites.add_row(position=56, ancestral_state='G')
tables.sites.add_row(position=64, ancestral_state='T')
tables.sites.add_row(position=72, ancestral_state='A')
tables.sites.add_row(position=80, ancestral_state='T')
tables.sites.add_row(position=88, ancestral_state='A')
tables.sites.add_row(position=96, ancestral_state='C')


### MUTATIONS TABLE ###
tables.mutations.add_row(site=2, node=0, derived_state='A')
tables.mutations.add_row(site=4, node=6, derived_state='C')
tables.mutations.add_row(site=6, node=5, derived_state='A')
tables.mutations.add_row(site=9, node=3, derived_state='A')
tables.mutations.add_row(site=11, node=0, derived_state='T')


#sort tables to make sure edges and sites tables are in the correct order
tables.sort()

fig2_ts = tables.tree_sequence()

fig2_ts.draw_svg(y_axis=True)


### NOTES ###
#assign example positions to breakpoint locations
#L0 0.0
#L1 30.0
#L2 60.0
#L3 100.0

#assign the node labels in figure (letters) to integer values:
#A 0
#B 1
#C 2
#D 3
#K 4
#P 5
#R 6
#W 7
#X 8

#assign example times to nodes:
#A T0 0
#B T0 0
#C T0 0
#D T0 0
#K TII 5
#P TIII 10
#R TV 15
#W TVII 20
#X TIX 25




#assign example positions for each site
#0 S0 8 G
#1 S1 16 T
#2 S2 24 G
#3 S3 32 G
#4 S4 40 T
#5 S5 48 C
#6 S6 56 G
#7 S7 64 T
#8 S8 72 A
#9 S9 80 T
#10 S10 88 A
#11 S11 96 C

#assign example times to each mutation
#0 2 0 A TVI 18
#1 4 6 C TVIII 22
#2 6 5 A TIV 13
#3 9 3 A TII 5
#4 11 0 T TI 3