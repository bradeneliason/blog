using Plots, DifferentialEquations, LinearAlgebra

## Creating a plotting recipe for these pursuit curves
@userplot PursuitPlot
@recipe function f(cp::PursuitPlot)
    t, u, i = cp.args
    u = hcat(u...)';
    x, y = u[1:i,1], u[1:i,2];
    # n = length(x)
    linewidth --> range(0.75, 10, length = i+1)
    aspect_ratio --> 1
    linealpha --> range(0.25, 1, length = i+1)
    x, y
end

## Define the pursuit differential equations
pursuit(u, k, t) = k * (pursued(t) - u) / norm(pursued(t)-u)

## Pursued Running in a Circle
pursued(t) = @. [cos(t), sin(t)]
prob = ODEProblem(pursuit, [4.0, 0.0], (0.0, 17), 0.8)
sol = solve(prob, saveat=0.1);
N = length(sol.t)

pursuitplot(sol.t, pursued.(sol.t), N, label="Rabbit")
pursuitplot!(sol.t, sol.u, N, label="Fox")
savefig("Images/PursuitCurves_Circle_R01.png")

# Generate Animation
anim = @animate for i ∈ 1:N
    pursuitplot(sol.t, pursued.(sol.t), i, label="Rabbit")
    pursuitplot!(sol.t, sol.u, i, label="Fox")
    plot!(xlims=(-1.1, 2.1), ylims=(-1.1,1.1))
end
gif(anim, "Images/PursuitCurves_Circle_R01.gif", fps = 15)

## Pursued Running in a Line
pursued(t) = @. [0, (t - 3)]
prob = ODEProblem(pursuit, [4.0, 0.0], (0.0,6), 1.0)
sol = solve(prob, saveat=0.1);
N = length(sol.t)

pursuitplot(sol.t, pursued.(sol.t),N, label="Rabbit")
pursuitplot!(sol.t, sol.u, N, label="Fox")
savefig("Images/PursuitCurves_Line_R01.png")

# Generate Animation
anim = @animate for i ∈ 1:N
    pursuitplot(sol.t, pursued.(sol.t), i, label="Rabbit")
    pursuitplot!(sol.t, sol.u, i, label="Fox")
    plot!(xlims=(-3,7), ylims=(-3.1,3.1))
end
gif(anim, "Images/PursuitCurves_Line_R01.gif", fps = 15)