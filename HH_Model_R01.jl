using Plots, DifferentialEquations, LaTeXStrings

## Constants
const C_m     = 1       # μF/cm², membrane capacitance
const G_NaMax = 120     # mS/cm², maximum conductivity of Na channel
const G_KMax  = 36      # mS/cm², maximum conductivity of K channel
const G_L     = 0.3     # mS/cm², leak conductivity
const V_r  = -65        # mV, resting potential
const V_Na =  50        # mV, Nernst voltage for Na 
const V_K  = -77        # mV, Nernst voltage for K
const V_L  = -54.387    # mV, Nernst voltage for leak

## Function Definitions
α_n(dV) = @. (0.1 - 0.01*dV)/(exp(1 - 0.1*dV)- 1)
β_n(dV) = @. 0.125/(exp(0.0125*dV))
α_m(dV) = @. (2.5 - 0.1*dV)/(exp(2.5 - 0.1*dV)- 1)
β_m(dV) = @. 4/(exp(dV/18))
α_h(dV) = @. 0.07/(exp(0.05*dV))
β_h(dV) = @. 1/(exp(3 - 0.1*dV)+ 1)

n_∞(dV) = @. α_n(dV)/(α_n(dV)+β_n(dV))
m_∞(dV) = @. α_m(dV)/(α_m(dV)+β_m(dV))
h_∞(dV) = @. α_h(dV)/(α_h(dV)+β_h(dV))

τ_n(dV) = 1 ./ (α_n(dV) .+ β_n(dV))
τ_m(dV) = 1 ./ (α_m(dV) .+ β_m(dV))
τ_h(dV) = 1 ./ (α_h(dV) .+ β_h(dV))

## Transfer Rate Coefficients Plots
p1 = plot([α_n β_n], -50, 50, lw=2, label=[L"\alpha_n" L"\beta_n"], legend=:top)
p2 = plot([α_m β_m], -50, 50, lw=2, label=[L"\alpha_m" L"\beta_m"])
p3 = plot([α_h β_h], -50, 50, lw=2, label=[L"\alpha_h" L"\beta_h"], legend=:top)
plot(p1, p2, p3, layout = grid(3, 1,), size=(500,500))
savefig("Images/HHModel_TransferRatePlot_R01.png")

p1 = plot([τ_n n_∞], -50, 50, lw=2, label=[L"\tau_n" L"n_{\inf}"])
p2 = plot([τ_m m_∞], -50, 50, lw=2, label=[L"\tau_m" L"m_{\inf}"], legend=:topleft)
p3 = plot([τ_h h_∞], -50, 50, lw=2, label=[L"\tau_h" L"h_{\inf}"])
plot(p1, p2, p3, layout = grid(3, 1,), size=(500,500))
savefig("Images/HHModel_SteadyState_TimeConstant_R01.png")

## Model Definition
I_inj(t) = (5 < t) & (t < 30) ? 10 : 0  # Injected Current Function

function HH_model(du,u,p,t)
    n, m, h, Vm = u 

    V_diff = Vm - V_r          # difference between the rest voltage and current membrane voltage
    G_K  = G_KMax  * n^4       # Sodium conductance
    G_Na = G_NaMax * h * m^3   # Potasium conductance
    
    # Update transfer rate coefficients, n, m, and h
    du[1] = α_n(V_diff)*(1-n) - β_n(V_diff)*n
    du[2] = α_m(V_diff)*(1-m) - β_m(V_diff)*m
    du[3] = α_h(V_diff)*(1-h) - β_h(V_diff)*h
    
    # Update cell membrane voltage, Vm
    du[4] = (I_inj(t) + (V_Na - Vm)*G_Na + (V_K - Vm)*G_K + (V_L - Vm)*G_L)/C_m   
end

## Run Model:
u0 = [n_∞(0); m_∞(0) ; h_∞(0); -65.1]
tspan = (0.0,50.0)
prob = ODEProblem(HH_model, u0, tspan)
sol = solve(prob, saveat=0.01);

## Plotting
p1 = plot(sol.t, sol[4,:], legend=false, lw=2, ylabel="Volatge [mV]")
p2 = plot(sol.t, I_inj.(sol.t), legend=false, lc=:red, lw=2, ylabel="Current")
p3 = plot(sol.t, sol[1:3,:]', label=["n" "m" "h"], legend=:topright, lw=2, 
        xlabel="Time [ms]", ylabel="Fraction Active")

l = grid(3, 1, heights=[0.4, 0.2 ,0.4])
plot(p1, p2, p3, layout = l, size=(800,500))
savefig("Images/HHModel_SimulationResult_R01.png")