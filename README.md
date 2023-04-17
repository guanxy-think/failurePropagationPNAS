# Code and data for the failure propagation study
## null model generation.R
This script generates the null models used by the inference model. Lines 7 to 49 gnerate the null models for the case study of congestion propagation in Sioux Falls nextwork. The generation process takes Sioux-Falls network topology as input (in the folder "benchmark model"). The generated null models preserve the degree of each node in Sioux-Falls network. Lines 60 to 155 generate the null models for the case study of failure propagations in interdependent transportation and power networks. The inputs are degree distributions of the power network, transportation networks and the interdependencies. The generated null models preserve the three degree distributions, respectievely.

## F-W.R
This script runs a traffic assignment model in the Sioux-Falls network to obtain the time when each road segment gets congested (i.e., the failure time). The inputs include the Sioux-Falls network topology and OD demand, both in the folder "benchmark model".
