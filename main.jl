## Main.jl
# MK(July 1, 2025, minkim.econ@gmail.com)

## Setting up the environment
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("replication_code.jl")

## Plots figures
plot_Ramsey_Markov(parameters) #fig1
plot_selfenforcingplan(parameters) #fig2
plot_sustainable(parameters)    #fig3
plot_compartivestatics(parameters) #fig4
plot_compartivestatics_path(parameters) #fig5
