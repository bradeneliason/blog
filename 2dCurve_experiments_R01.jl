# using ApproxFun
using DifferentialEquations
using LinearAlgebra
using Plots

## Swimming Dog curves
# https://mathcurve.com/courbes2d.gb/nageur/nageur.shtml
# ẋ = -Vn * x /√x²+y²
# ẏ = Vc - Vn * y /√x²+y²

# function swimmingdog(u, p, t)
#     x, y = u
#     dx = -x /norm(u)
#     dy = p - y /norm(u)
#     [dx, dy]
# end

swimmingdog(u, p, t) = [-u[1]/norm(u); p-u[2]/norm(u)]

u0 = [1.0; 0.0]
tspan = (0.0,10.0)

sd_prob_06 = ODEProblem(swimmingdog, u0, tspan, 0.6)
sd_prob_08 = ODEProblem(swimmingdog, u0, tspan, 0.8)
sd_prob_10 = ODEProblem(swimmingdog, u0, tspan, 1.0)
sd_prob_12 = ODEProblem(swimmingdog, u0, tspan, 1.2)
sd_prob_14 = ODEProblem(swimmingdog, u0, tspan, 1.4)

sd_sol_06 = solve(sd_prob_06, saveat=0.1);
sd_sol_08 = solve(sd_prob_08, saveat=0.1);
sd_sol_10 = solve(sd_prob_10, saveat=0.1);
sd_sol_12 = solve(sd_prob_12, saveat=0.1);
sd_sol_14 = solve(sd_prob_14, saveat=0.1);

# Plotting
plot( sd_sol_06, vars=(1,2), label="0.6", frame=:zerolines, lw=3, lc=:red)
plot!(sd_sol_08, vars=(1,2), label="0.8", frame=:zerolines, lw=3, lc=:red)
plot!(sd_sol_10, vars=(1,2), label="1.0", frame=:zerolines, lw=3, lc=:green)
plot!(sd_sol_12, vars=(1,2), label="1.2", frame=:zerolines, lw=3, lc=:blue)
plot!(sd_sol_14, vars=(1,2), label="1.4", frame=:zerolines, lw=3, lc=:blue)
plot!(ylims=(0,2.0), aspect_ratio=:equal, 
    title="Swimming Dog Curves", size = (500,800)
)

## Tractrix
# https://mathcurve.com/courbes2d.gb/tractrice/tractrice.shtml
# function tractrix(u, p, t)
#     x, y = u
#     dx = 1
#     dy = -y /norm([p,y])
#     [dx, dy]
# end

tractrix(u, p, t) = [1, -u[2] /norm([p,u[2]])]

u0 = [0.0; 1.0]
tspan = (0.0,30.0)

tx_prob_01 = ODEProblem(tractrix, u0, tspan, 1.0)
tx_prob_02 = ODEProblem(tractrix, u0, tspan, 2.0)
tx_prob_05 = ODEProblem(tractrix, u0, tspan, 5.0)
tx_prob_10 = ODEProblem(tractrix, u0, tspan, 10.0)
tx_prob_20 = ODEProblem(tractrix, u0, tspan, 20.0)

tx_sol_01 = solve(tx_prob_01, saveat=0.1);
tx_sol_02 = solve(tx_prob_02, saveat=0.1);
tx_sol_05 = solve(tx_prob_05, saveat=0.1);
tx_sol_10 = solve(tx_prob_10, saveat=0.1);
tx_sol_20 = solve(tx_prob_20, saveat=0.1);

# Plotting
plot( tx_sol_01, vars=(1,2), label="a = 1", frame=:zerolines, lw=3)
plot!(tx_sol_02, vars=(1,2), label="a = 2", frame=:zerolines, lw=3)
plot!(tx_sol_05, vars=(1,2), label="a = 5", frame=:zerolines, lw=3)
plot!(tx_sol_10, vars=(1,2), label="a = 10", frame=:zerolines, lw=3)
plot!(tx_sol_20, vars=(1,2), label="a = 20", frame=:zerolines, lw=3)
plot!(ylims=(0,1.2), title="Tractrix")


## Elastic catenary  -- TODO
# https://mathcurve.com/courbes2d.gb/chainette/chainetteelastique.shtml
catenary(dy, y, k, t) = ddy = sqrt(1 + dy^2) / (1 + k*sqrt(1+dy^2))

u0 = 0.0
du0= -1.0
tspan = (0.0, 2.0)

ct_prob_01 = SecondOrderODEProblem(catenary, du0, u0, tspan, 0.2)
ct_sol_01 = solve(ct_prob_01, saveat=0.01);

plot(ct_sol_01, vars=(1,2), label="a = 1", frame=:zerolines, lw=3)



## Elastica -- TODO
# https://mathcurve.com/courbes2d.gb/linteaire/linteaire.shtml
# function elastica(du, u, a, t)
#     x, y = u
#     dx, dy = du
#     dx = 1
#     ddy = 2x * (1+ dy^2)^1.5 / a
#     [1.0; ddy]
# end
function elastica(u, p, t)
    x, y = u
    k, a = p
    dx = 1
    # dy =(x^2 - k*a^2)/sqrt(a^4 - (x^2-k*a^2)^2)
    dy = (x^2/a^2) - k
    [1.0; dy]
end

# elastica(u, p, t) = [1, -u[2] /norm([p,u[2]])]

u0 = [0.0; 1.0]
p = [0.2, 1.5]
du0= [0.0; -1.0]
tspan = (0.0, 10.0)

el_prob_01 = SecondOrderODEProblem(elastica, du0, u0, tspan, 1.0)
# el_prob_01 = ODEProblem(elastica, u0, tspan, p)
el_sol_01 = solve(el_prob_01, saveat=0.01);
el_sol_01




## DELAUNAY ROULETTE
# https://mathcurve.com/courbes2d.gb/delaunay/delaunay.shtml
function roulette(dy, y, p, t)
    a, b = p
    ε = b/a
    tmp = (y^2 + ε*b^2)
    ddy = (4*a^2*y^2 - tmp^2) / tmp^2
end

u0 = 1.0
du0= 1.0
p = [1.0, 1.2]
tspan = (0.0, 10.0)

dr_prob_01 = SecondOrderODEProblem(roulette, du0, u0, tspan, p)
dr_sol_01 = solve(dr_prob_01,);

plot(dr_sol_01, vars=(1,2), 
    frame=:zerolines,
    lw=3,
    legend=:none
)
