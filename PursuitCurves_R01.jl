using DifferentialEquations
using LinearAlgebra
using Plots

## Creating a plotting recipe for these pursuit curves
@userplot PursuitPlot
@recipe function f(cp::PursuitPlot)
    t, u, i = cp.args
    u = hcat(u...)';
    x, y = u[1:i,1], u[1:i,2];
    linewidth --> range(0.75, 10, length = i+1)
    aspect_ratio --> 1
    linealpha --> range(0.25, 1, length = i+1)
    x, y
end

plotemoji!(u, symbol) = annotate!(u[1], u[2], Plots.text(symbol, "courier"))

# Define the pursuit differential equations
pursuit(u, k, t) = k * (pursued(t) - u) / norm(pursued(t)-u)

## Pursued Running in a Circle
pursued(t) = [cos(t), sin(t)]
prob = ODEProblem(pursuit, [4.0, 0.0], (0.0, 17), 0.8)
sol = solve(prob, saveat=0.1);
N = length(sol.t)

pursuitplot(sol.t, pursued.(sol.t), N, label="Rabbit")
# plotemoji!(pursued(sol.t[end]), "🐰")
pursuitplot!(sol.t, sol.u, N, label="Fox")
# plotemoji!(sol.u[end], "🦊")
savefig("./images/PursuitCurves_Circle_R02.png")
# Error with saving plot as a png
# GKS: glyph missing from current font

## Generate Animation
anim = @animate for i ∈ 1:N
    pursuitplot(sol.t, pursued.(sol.t), i, label="Rabbit")
    pursuitplot!(sol.t, sol.u, i, label="Fox")
    plot!(xlims=(-1.1, 2.1), ylims=(-1.1,1.1))
end
gif(anim, "./images/PursuitCurves_Circle_R02.gif", fps = 15)

## Pursued Running in a Line
pursued(t) = [0, (t - 3)]
prob = ODEProblem(pursuit, [4.0, 0.0], (0.0,6), 1.0)
sol = solve(prob, saveat=0.1);
N = length(sol.t)

pursuitplot(sol.t, pursued.(sol.t),N, label="Rabbit")
pursuitplot!(sol.t, sol.u, N, label="Fox")
savefig("./images/PursuitCurves_Line_R02.png")

## Generate Animation
anim = @animate for i ∈ 1:N
    pursuitplot(sol.t, pursued.(sol.t), i, label="Rabbit")
    pursuitplot!(sol.t, sol.u, i, label="Fox")
    plot!(xlims=(-3,7), ylims=(-3.1,3.1))
end
gif(anim, "./images/PursuitCurves_Line_R02.gif", fps = 15)