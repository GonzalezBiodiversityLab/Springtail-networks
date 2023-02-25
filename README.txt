Description of variables used in statistical analysis


# General

Config - habitat configuration (1 = lattice, 2 = partially random, 3 = random); note that partially random networks are referred to as 'small world' (SMW) throughout the script.
Rep - replicate network within each configuration; 8 replicates of each configuration (Rep 2 - 9) - note: for lattice configuration in the model, there is only one Rep (Rep = 2) because there is only one topology possible and no physical "replicates". (Note: Reps A-H in Figure S2 correspond to Reps 2-9 in data files).
node - indicates node within network; each network has ten nodes, 1 - 10.
landscapeID (experiment only) - unique physical network (unique Rep/Config combination).
topology (model only) - network topology (unique Rep/Config combination).
ModelRep - (model only) simulations are run multiple times on the same network topologies; ModelRep indicates the model simulation run.
repID (model only) = unique simulated network (unique ModelRep/Rep/Config combination).
Day - days since start of experiment; Ranges from 1-182 in experiment and 1 - 200 in model.


# Spatial properties of nodes and networks

distToSource - minimum distance in mm through links from source node (node 5 in all networks) to focal node.
AlgebraicConnectivity_Binary - algebraic connectivity of the network
AlgebraicConnectivity_Weighted - algebraic connectivity of the network, weighted by dispersal probability (a function of link length).
NetworkDiameter_Binary - network diameter in links
NetworkDiameter_Weighted - network diameter in mm

# Experiment and model population size data

N - in 'long' data table format, gives observed number of collembola in a node.
Node1, Node2, Node3, etc. - in 'wide' data table format, gives observed number of collembola in a node.

# General model

B = shape of the dispersal kernel; exponent of the negative exponential function fit to the dispersal kernel (dispersal probability as a function of link length).
Clnz.Time = number of days to node colonization (first day on which node population size is observed to be 1 or greater).
Occp.Time = number of days to full network occupancy (first day on which observed population size of all nodes in a network is 1 or greater).