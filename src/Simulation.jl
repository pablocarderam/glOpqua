using DifferentialEquations

function simulate(
        params::ModelParameters, init::Vector{Float64}, tspan::Tuple{Float64, Float64};
        odeFunc::Function=glOpqua.multistrain2Immunity, solver=QNDF(autodiff=false)
        )
    ode_prob = ODEProblem(odeFunc, init, tspan, params)
    return solve(ode_prob, solver, saveat=0.1)
end
