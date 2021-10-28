# using ApproxFun
using DifferentialEquations
using LinearAlgebra
using Plots

function solvemany(func, u0, tspan, params, saveat=0.1)
    solns = []
    for p in params
        prob = ODEProblem(func, u0, tspan, p)
        push!(solns, solve(prob, saveat=saveat))
    end
    return solns
end

## Swimming Dog curves
# https://mathcurve.com/courbes2d.gb/nageur/nageur.shtml
swimmingdog(u, p, t) = [-u[1]/norm(u); p-u[2]/norm(u)]

u0 = [1.0; 0.0]
tspan = (0.0,10.0)
params =  [0.6, 0.8, 1.0, 1.2, 1.4]
sd_solns = solvemany(swimmingdog, u0, tspan, params);

plot()
cols = [:red, :red, :green, :blue, :blue]
for (s, p, c) in zip(sd_solns, params, cols)
    plot!(s, vars=(1,2), label="$p", lw=3, lc=c)
end
plot!(title="Swimming Dog Curves", 
    ylims=(0,1.0), 
    frame=:zerolines,
    aspect_ratio=:equal,
    size=(500,500),
)
savefig("./images/swimmingdog_curves_R01.png")

## Tractrix
# https://mathcurve.com/courbes2d.gb/tractrice/tractrice.shtml
tractrix(u, p, t) = [1, -u[2] /norm([p,u[2]])]

u0 = [0.0; 1.0]
tspan = (0.0,30.0)
params =  [1.0, 2.0, 5.0, 10.0, 20.0]
tx_solns = solvemany(tractrix, u0, tspan, params);

plot()
for (s, p) in zip(tx_solns, params)
    plot!(s, vars=(1,2), label="$p", lw=3)
end
plot!(title="Tractrix Curves", xlims=(0,30), frame=:zerolines)
savefig("./images/tractrix_curves_R01.png")
