## Solving for a Uniformly loaded beam in 7 lines of Julia
using ApproxFun, Plots
L, E, I = 12.0, 1.0, 1.0
d = 0..L
z = Fun(identity, d)
B = [[Evaluation(d,leftendpoint,k) for k=0:1]... ; [Evaluation(d,rightendpoint,k) for k=2:3]... ;]
v = [B; E*I*Derivative()^4] \ [ zeros(4)..., -one(z)]
func_name = zip([v, v', v'', v'''], ["Deflection", "Angle", "Momement", "Shear"])
plot([plot(z, f, title=n, label="") for (f,n) in func_name]..., lw=3)


## Okay Solving it the normal way
using ApproxFun, Plots
° = π/180;

# Setting up problem
L = 2.0                 # Length in m
E = 82.74e9             # Elasticin Pa
I = 444.0 * (0.01)^4    # Moment of interia m⁴
d = 0..L                # Domain of beam
z = Fun(identity, d)    # Length along beam in m
D = Derivative()

# Problem Definition
# q = -1000*sin(2π*z/L)   # Sinusoidally loaded in N/m
# q = 1000*(z/L)          # Triangular loading in N/m
q = 1000*one(z)         # Uniformly loading in N/m
w = E*I*D^4             # DiffEq for beam deflection

# Boundary Conditions
B= [Evaluation(d, 0, 0) # Beam vertically constrained at z = 0
    Evaluation(d, 0, 1) # Beam's angle is constrained at z = 0
    Evaluation(d, L, 2) # Beam's moment is 0  at z = L
    Evaluation(d, L, 3)]# Beam's shear is 0 at z = L

# Solving for vertical displacement
v = [B; w] \ [ zeros(4)...; -q]

# Renaming and scaling variables
θ, M, V = (v'/°),(v''*E*I), (v'''*E*I)

println("Max Deflection: \t $(round(1000*v(L),digits=2)) mm")
println("Max Angle:      \t $(round(θ(L),digits=3))° ")

# Plotting
p1 = plot(z, 1000v, legend=:none, title="Deflection [mm]", lc=:blue)
p2 = plot(z, θ, legend=:none, title="Angle [°]", lc=:orange)

plot(z, M/1000, label="Momement [kN⋅m]",
    fill = (0, 0.15, :blue),  lc=:blue)
plot!(z, V/1000, label="Shear [kN]",
    fill = (0, 0.15, :green), lc=:green)
p3 = plot!(z, -q/1000, label="Load [kN/m]", # line=:stem,
    fill = (0, 0.15, :red),   lc=:red)

l = grid(3, 1, heights=[0.2, 0.2 ,0.6])
plot(p1, p2, p3, layout=l,  lw = 3, size=(500, 900), frame=:zerolines)
# savefig("Beam_deflection_plot_R01.png")