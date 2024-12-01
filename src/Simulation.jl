using DifferentialEquations
using LSODA

function simulate(
    params::ModelParameters, init::Vector{Float64}, tspan::Tuple{Float64,Float64};
    odeFunc::Function=glOpqua.multistrainImprintMature, solver="None", saveat=0.1
)
    ode_prob = ODEProblem(odeFunc, init, tspan, params)

    if solver == "None"
        return solve(ode_prob, saveat=0.1)
    elseif solver == "LSODA"
        return solve(ode_prob, lsoda(), saveat=saveat)
    elseif solver == "QNDF"
        return solve(ode_prob, QNDF(autodiff=false), saveat=0.1)
    else
        return solve(ode_prob, solver, saveat=0.1)
    end
end
