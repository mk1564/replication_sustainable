## Sustainable FXI
# June 29, 2025
# Min Kim (minkim.econ@gmail.com)

using Parameters
using Plots, LaTeXStrings #, QuantEcon, LinearAlgebra
# import Random
# Random.seed!(1234)

## Setting up parameter values
parameters = @with_kw (α=0.8, # responsiveness of capital flow to exchange rate
                        θ=5, # cost of FX reserve management
                        γ=0.15, # responsiveness of current account to exchange rate
                        ebar=0.5, # bliss point of exchange rate
                        T=1000, Ts=4, fsbar=0.01, T_Plot=41, 
                        θs = [1, 3, 5, 7, 10])
                        #σ=0.5)


function solve_Ramsey_Markov(parameters)
    @unpack α, θ,  γ, ebar, T, Ts, fsbar, T_Plot, θs = parameters() 
    β = α/(α+ γ) # discount factor
    U(e,f) = ((e .- ebar).^2/2 .+ θ/2*f.^2) # one-period welfare

    ## Ramsey plan
    κ = 1/(α+ γ)
    θ_inv = 1/θ 
    λ = θ_inv * κ^2
    ω_R = (1 + β + λ - sqrt((1+β+λ)^2 - 4*β))/(2*β)
    δ = ebar * (1-β) / κ
    e0_R = (β*(1-ω_R) + λ)/(1-β*ω_R + λ) * ebar # 

    # Compute time path for Ramsey plan
    e_Ramsey = zeros(T)
    e_Ramsey[1] = e0_R # initial value
    for t in 2:T
        e_Ramsey[t] = ω_R * e_Ramsey[t-1]  + (1 - ω_R) * ebar 
    end
    f_Ramsey = (α+ γ) * (1 - ω_R * β) .* (e_Ramsey .- ebar) .+ γ * ebar
    Jhat_Ramsey_function = x -> - (1-β) * ebar^2 / κ^2 .- (2*(1-β)*ebar)/κ^2 .* x .-  ((1-β*ω_R)/κ^2 + 1/θ) .* x.^2
    # write a function to compute Jhat_Ramsey using x as input
    J_Ramsey =  θ/2 * Jhat_Ramsey_function(e_Ramsey .- ebar) # compute Jhat_Ramsey for each e_Ramsey
    # J_Ramsey2 = -U(e_Ramsey,f_Ramsey)/(1-β) # compute J_Ramsey
    # Compute time path for Markov
    ω_MPE = 1/(1+θ* γ*(α+ γ)) 
    e_MPE = ω_MPE * ebar * ones(T)
    f_MPE = γ * ω_MPE * ebar * ones(T)
    J_MPE = -U(e_MPE,f_MPE)/(1-β)

    # Equilibrium path for Ramsey and Markov
    f = [f_Ramsey f_MPE] 
    e = [e_Ramsey e_MPE] 
    J = [J_Ramsey J_MPE] 
    Δk = f .- γ .* e # change in capital
    ca = γ .* e # current account

    return f, e, J, Δk, ca
end

function plot_Ramsey_Markov(parameters)
    @unpack α, θ,  γ, ebar, T, Ts, fsbar, T_Plot, θs = parameters() 
    f, e, J, Δk, ca = solve_Ramsey_Markov(parameters)

    p1 = plot(e[1:T_Plot,:], label=["Ramsey" "Markov"], xlabel="Time", ylabel="", title=L"e_t", legend=:topright, fg_legend = :false)
    p2 = plot(f[1:T_Plot,:], xlabel="Time", ylabel="", title=L"f_t", legend=:false, fg_legend = :false)
    p3 = plot(Δk[1:T_Plot,:], xlabel="Time", ylabel="", title=L"\Delta k_t", legend=:false, fg_legend = :false)
    p4 = plot(J[1:T_Plot,:], xlabel="Time", ylabel="", title=L"Welfare", legend=:false, fg_legend = :false)
    fig = plot(p1,p2,p3,p4)
    savefig(fig,"fig1.pdf")
    # display(fig)
end

# plot_Ramsey_Markov(parameters)

function solve_sustainableplan(parameters)
    @unpack α, θ,  γ, ebar, T, Ts, fsbar, T_Plot, θs = parameters() 
    β = α/(α+ γ) # discount factor
    U(e,f) = ((e .- ebar).^2/2 .+ θ/2*f.^2) # one-period welfare

    f, e, J, Δk, ca = solve_Ramsey_Markov(parameters)
    f_Ramsey = f[:,1] # Ramsey f
    e_Ramsey = e[:,1] # Ramsey e

    fs = [ones(Ts)*fsbar; f_Ramsey]
    es = zeros(Ts+T)

    # calculate discount sum
    discount = ones(Ts+T)
    for j in 2:Ts+T # j=1 is 1
        discount[j] = (α/(α+ γ))^(j-1) #* f_series[j]
    end

    for t in 1:Ts+T
        es[t] = 1/(α+ γ) * sum(discount[1:end-(t-1)] .*fs[t:end]) # forward-looking 
    end

    # Calculate utility of stick plan
    Us = zeros(Ts+T) # in each period
    for t in 1:Ts+T
        Us[t] = β^(t-1) * -U(es[t],fs[t])
    end

    Vs = zeros(Ts+T) # value
    for t in 1:Ts+T
        Vs[t] = sum(Us[t:end] / β^(t-1)) # divide by β^(t-1) to make discount factor to be 1 at each period t
    end

    # value from deviation 
    Vd = zeros(Ts+T)
    for t in 1:Ts+T
        Vd[t] = -U(es[t],0) + β * Vs[1]
    end
    
    # value by devating from sustainable plan
    J_deviate = zeros(1,T)
    for t in 1:T
        J_deviate[t] = -U(e_Ramsey[t],0)  + β * Vs[1]
    end

    return es, fs, Vs, Vd, J_deviate
end

function plot_selfenforcingplan(parameters)
    @unpack α, θ,  γ, ebar, T, Ts, fsbar, T_Plot, θs = parameters()
    es, fs, Vs, Vd, J_deviate = solve_sustainableplan(parameters)
    l = @layout [α  θ;  γ]
    p1 = plot([0:T_Plot-1],es[1:T_Plot], title = L"e_t^S", label = "", legend= :bottomright)
    p2 = plot([0:T_Plot-1],fs[1:T_Plot], title = L"f_t^S", label = "", legend= :bottomright)
    p3 = plot([0:T_Plot-1],Vs[1:T_Plot], label = "Self-enforcing plan", title = L"Welfare", 
    #legend=(0.8,0.6), 
    legend = :bottomright,fg_legend = :false)
    plot!(p3,[0:T_Plot-1],Vd[1:T_Plot], label = "Deviation")
    # display(plot(p1,p2))
    # display(p3)
    p_combined = plot(p1,p2,p3, 
        layout = l, xlabel = "time")
    savefig(p_combined,"fig2.pdf")
    # display(p_combined)

    # check if higher welfare
    T_Plot == sum(Vs[1:T_Plot] .>  Vd[1:T_Plot])
end

# plot_selfenforcingplan(parameters)

function plot_sustainable(parameters)
    @unpack α, θ,  γ, ebar, T, Ts, fsbar, T_Plot, θs = parameters()
    
    f, e, J, Δk, ca = solve_Ramsey_Markov(parameters) 
    J_Ramsey = J[:,1] # Ramsey J
    es, fs, Vs, Vd, J_deviate = solve_sustainableplan(parameters)
    plt = plot([0:T_Plot-1],J_Ramsey[1:T_Plot], 
        xlabel = "time", title = L"Welfare", label="Sustainable plan", fg_legend = :false)
    plot!(plt, [0:T_Plot-1],J_deviate'[1:T_Plot], label="Deviation")
    # display(plt)
    savefig(plt,"fig3.pdf")
end
# plot_sustainable(parameters)

## Plots in Appendix 
function plot_compartivestatics(parameters)
    @unpack α, θ,  γ, ebar, T, Ts, fsbar, T_Plot, θs = parameters() 
    β = α/(α+ γ) # discount factor
    U(e,f) = ((e .- ebar).^2/2 .+ θ/2*f.^2) # one-period welfare
    κ = 1/(α+ γ)
    λ = 1/θ  * κ^2
    ω = (1 + β + λ - sqrt((1+β+λ)^2 - 4*β))/(2*β)
    eR0 = (β*(1-ω) + λ)/(1-β*ω + λ) * ebar
    e_grid = range(eR0-0.05, ebar + ebar - eR0, length = 100)
    ep_function =  ω.*e_grid .+ (1-ω)* ebar
    # looping over different θ
    labels = ["θ = $θ" for θ in θs]
    ep_functions_θ = zeros(length(θs),100)
    f_functions_θ = zeros(length(θs),100)
    eR0_functions_θ = zeros(length(θs))
    ω_functions_θ = zeros(length(θs))
    for (i,θ) in enumerate(θs)
        λ = 1/θ  * κ^2
        ω = (1 + β + λ - sqrt((1+β+λ)^2 - 4*β))/(2*β)
        eR0 = (β*(1-ω) + λ)/(1-β*ω + λ) * ebar
        ep_function =  ω.*e_grid .+ (1-ω)* ebar
        ω_functions_θ[i] = ω
        f_function = ((α+ γ)*(1-ω*β)*(e_grid .- ebar) .+ γ*ebar)
        f_functions_θ[i,:] = f_function
        ep_functions_θ[i,:] = ep_function
        eR0_functions_θ[i] = eR0
    end

    fig_policy_ep_theta = plot(e_grid, ep_functions_θ', label = hcat(labels...), xlabel = L"e", title = L"e'", fg_legend = :false, xlim = (e_grid[1],0.6), legend = :topleft)
    plot!(e_grid, e_grid, label="45 degree line", fg_legend = :false, xlabel = L"e", linestyle = :dash, color = :black, alpha = 0.5)

    fig_policy_f_theta = plot(e_grid, f_functions_θ', label = hcat(labels...), xlabel = L"e", title = L"f", fg_legend = :false, xlim = (e_grid[1],0.6))
    # scatter!([eR0_functions_θ], [((α+ γ).*(1 .- ω_functions_θ.*β).*(eR0_functions_θ .- ebar) .+ γ*ebar)], label = "", fg_legend = :false, marker = :circle, color = :black, alpha = 0.5)
    vline!([ebar], label = "", fg_legend = :false, linestyle = :dash, color = :black, alpha = 0.5)

    fig_policy_theta = plot(fig_policy_ep_theta, fig_policy_f_theta)
    savefig(fig_policy_theta,"fig4.pdf")
    # display(fig_policy_theta)
end

# plot_compartivestatics(parameters)

function plot_compartivestatics_path(parameters)
    @unpack α, θ,  γ, ebar, T, Ts, fsbar, T_Plot, θs = parameters() 
    β = α/(α+ γ) # discount factor
    U(e,f) = ((e .- ebar).^2/2 .+ θ/2*f.^2) # one-period welfare

    κ = 1/(α+ γ)
    λ = 1/θ  * κ^2
    # loop over different θ
    labels = ["θ = $θ" for θ in θs]
    T_Plot2 = 20 # shorter time horizon for plotting
    e_paths = zeros(length(θs),T_Plot2)
    f_paths = zeros(length(θs),T_Plot2)
    J_paths = zeros(length(θs),T_Plot2)
    u_paths = zeros(length(θs),T_Plot2)
    u_paths_fpart = zeros(length(θs),T_Plot2)
    u_paths_epart = zeros(length(θs),T_Plot2)
    for (i,θ) in enumerate(θs)
        λ = 1/θ  * κ^2
        ω = (1 + β + λ - sqrt((1+β+λ)^2 - 4*β))/(2*β)
        eR0 = (β*(1-ω) + λ)/(1-β*ω + λ) * ebar
        e_path = zeros(T_Plot2)
        e_path[1] = eR0
        for t in 2:T_Plot2
            e_path[t] = ω*e_path[t-1] + (1-ω)* ebar
        end
        f_path = ((α+ γ)*(1-ω*β)*(e_path .- ebar) .+ γ*ebar)
        J_path = θ/2*(-(1-β)*ebar^2/ κ^2 + 2*(1-β)*ebar/ κ^2 * ebar - ((1-β*ω)/κ^2 + 1/θ)*ebar^2) .-θ/2*(2*(1-β)*ebar/ κ^2) .+θ/2*(2*ebar*((1-β*ω)/κ^2 + 1/θ)).*e_path .- θ/2*((1-β*ω)/κ^2 + 1/θ).*e_path.^2
        u_path_fpart = -θ/2*f_path.^2
        u_path_epart = -(e_path .- ebar).^2/2
        
        
        
        e_paths[i,:] = e_path
        f_paths[i,:] = f_path
        J_paths[i,:] = J_path
        u_paths_fpart[i,:] = u_path_fpart
        u_paths_epart[i,:] = u_path_epart
        u_paths[i,:] = u_path_fpart + u_path_epart
    end     
        p_e_cs = plot(e_paths', xticks = 0:5:T_Plot2, legend=:bottomright, label = hcat(labels...), xlabel = "time", title = L"e_t", fg_legend = :false)
        
        p_f_cs = plot(f_paths', xticks = 0:5:2T_Plot2, legend=:bottomright, label = hcat(labels...), xlabel = "time", title = L"f_t", fg_legend = :false)
        
        p_cs = plot(p_e_cs ,p_f_cs)
        savefig(p_cs,"fig5.pdf")
        # display(p_cs)
end
# plot_compartivestatics_path(parameters)



# # κ * δ / (1- β*ω + λ)
# ω_R * ebar
# -F[1,1]

# 1/α*(α + γ + F[1,2]) # = omega

# 1/α*F[1,1] # = (1-ω_R) * ebar 

# -1/κ * β*(1-ω_R) * ebar # -F[1,1]
# 1/κ * (1-β*(ω_R)) # -F[1,2]

# (1-ω_R) * ebar 
# κ * δ / (1- β*ω_R + θ_inv*κ^2)


# -P[1,1]
# θ/2*(-(1-β)*ebar^2/ κ^2 + 2*(1-β)*ebar/ κ^2 * ebar - ((1-β*ω_R)/κ^2 + 1/θ)*ebar^2)
# # this is -P[1,1]

# -2*P[1,2]
# -θ/2*(2*(1-β)*ebar/ κ^2) +θ/2*(2*ebar*((1-β*ω_R)/κ^2 + 1/θ))
# # this is -2*P[1,2]

# -P[2,2]
# -θ/2*((1-β*ω_R)/κ^2 + 1/θ) # this is -P[2,2/]

# # compute this
# θ/2*(-(1-β)*ebar^2/ κ^2 + 2*(1-β)*ebar/ κ^2 * ebar - ((1-β*ω_R)/κ^2 + 1/θ)*ebar^2) + -θ/2*(2*(1-β)*ebar/ κ^2) +θ/2*(2*ebar*((1-β*ω_R)/κ^2 + 1/θ))
# + -θ/2*((1-β*ω_R)/κ^2 + 1/θ)  



# eR0
# -(1-β)*ebar/ ((1-β*ω_R) + λ) + ebar # this is eR0
# (β*(1-ω_R) + λ)/(1-β*ω_R + λ) * ebar # this is eR0

# ## Compute values
# function solve_Ramsey_Markov(T)
#     R = [ebar^2/2  -ebar/2;
#         -ebar/2 1/2]
#     Q = θ / 2
#     A = [1 0; 0 (α+ γ)/α]
#     B = [0; -1/α]
#     C = [1 0; 0  γ; 0 -γ]
#     D = [0;0;1]

#     ## Solve Ramsey 
#     #LQ Problem (Subproblem 1)
#     lq = QuantEcon.LQ(Q, R, A, B, bet=β)
#     P, F, d = stationary_values(lq)
#     # Solve Subproblem 2
#     eR0 = -P[1, 2] / P[2, 2]
#     x0 = [1; eR0]

#     ## Solve Markov
#     fMPE = ( γ/(1+θ* γ*(α+ γ)))*ebar
#     eMPE = (1/(1+θ* γ*(α+ γ)))*ebar
#     JMPE = -U(eMPE,fMPE)/(1-β)

#     ## Time domain
#     # Ramsey
#     x_sequence, f_sequence = compute_sequence(lq, x0, T)
#     x_sequence = x_sequence[:,1:end-1]

#     y_sequence = zeros(3,T)
#     J_sequence = zeros(1,T)
#     for t in 1:T
#     y_sequence[:,t] = C * x_sequence[:,t] + D .* f_sequence[:,t]
#     J_sequence[t] = -transpose(x_sequence[:,t]) * P * x_sequence[:,t]
#     end

#     e_Ramsey = x_sequence[2,:]
#     f_Ramsey = f_sequence'
#     ca_Ramsey = y_sequence[2,:]
#     Δk_Ramsey = y_sequence[3,:]
#     k_Ramsey = cumsum(Δk_Ramsey)
#     R_Ramsey = cumsum(f_sequence,dims=2)'
#     J_Ramsey = J_sequence'

#     e_MPE = eMPE .* ones(T,1)
#     f_MPE = fMPE .* ones(T,1)
#     ca_MPE =  γ .* e_MPE
#     Δk_MPE = ca_MPE .- f_MPE
#     k_MPE = cumsum(Δk_MPE, dims=1)
#     R_MPE = cumsum(f_MPE, dims=1)
#     J_MPE = JMPE .* ones(T,1)

#     return e_Ramsey, f_Ramsey, ca_Ramsey, Δk_Ramsey, k_Ramsey, R_Ramsey, J_Ramsey, e_MPE,  f_MPE, ca_MPE,     Δk_MPE, k_MPE, R_MPE, J_MPE
# end

# e_Ramsey_QE, f_Ramsey_QE, ca_Ramsey_QE, Δk_Ramsey_QE, k_Ramsey_QE, R_Ramsey_QE, J_Ramsey_QE, e_MPE_QE,  f_MPE_QE, ca_MPE_QE, Δk_MPE_QE, k_MPE_QE, R_MPE_QE, J_MPE_QE = solve_Ramsey_Markov(T) 

# ## Plots

# function plot_Ramsey_Markov(T_Plot)
  
#     p1=plot([0:T_Plot-1],e_Ramsey[1:T_Plot], 
#     xlabel = "time", title = L"e_t", label="Ramsey", fg_legend = :false)
#     plot!(p1, [0:T_Plot-1],e_MPE[1:T_Plot], label="Markov")
#     p2=plot([0:T_Plot-1],ca_Ramsey[1:T_Plot], 
#     xlabel = "time", title = L"ca_t", label="Ramsey", legend=:false)
#     plot!(p2, [0:T_Plot-1],ca_MPE[1:T_Plot], label="Markov")
#     p3=plot([0:T_Plot-1],Δk_Ramsey[1:T_Plot],
#     xlabel = "time",  title=L"\Delta k_t", label="Ramsey", legend=:false)
#     plot!(p3, [0:T_Plot-1],Δk_MPE[1:T_Plot], label="Markov")
#     p4=plot([0:T_Plot-1],f_Ramsey[1:T_Plot], 
#     xlabel = "time", title = L"f_t", label="Ramsey", legend=:false)
#     plot!(p4, [0:T_Plot-1],f_MPE[1:T_Plot], label="Markov")
#     Ramsey_Markov_plot1 = plot(p1,p2,p3,p4)
#     savefig(Ramsey_Markov_plot1,"Ramsey_Markov_plot1.pdf")
#     display(Ramsey_Markov_plot1)

#     ## capital and reserve (starting 0 initial values)
#     p5 = plot([0:T_Plot-1],k_Ramsey[1:T_Plot], 
#         xlabel = "time", title = L"k_t", label="Ramsey", fg_legend = :false)
#     plot!(p5,[0:T_Plot-1], k_MPE[1:T_Plot], label="Markov")
#     p6 = plot([0:T_Plot-1],R_Ramsey[1:T_Plot],
#         xlabel = "time", title = L"R_t", label="Ramsey", legend=:false)
#     plot!(p6, [0:T_Plot-1],R_MPE[1:T_Plot], label="Markov")
#     Ramsey_Markov_plot2 = plot(p5,p6)
#     savefig(Ramsey_Markov_plot2,"Ramsey_Markov_plot2.pdf")
#     display(Ramsey_Markov_plot2)

#     ## welfare
#     p7 = plot([0:T_Plot-1],J_Ramsey[1:T_Plot], 
#         xlabel = "time", title = "Welfare", label="Ramsey", fg_legend = :false)
#     plot!(p7, [0:T_Plot-1],J_MPE[1:T_Plot], label="Markov")
#     savefig(p7,"Ramsey_Markov_plot3.pdf")
#     display(p7)
# end

# plot_Ramsey_Markov(T_Plot)

# ## Self-enforcing plan and Sustainable plan

