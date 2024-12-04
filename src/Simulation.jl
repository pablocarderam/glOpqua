using DifferentialEquations
using LSODA
using Random
using StatsBase

function simulateODE(
    parameters::ModelParameters, init::Vector{Float64}, tspan::Tuple{Float64,Float64};
    odeFunc::Function=glOpqua.multistrainImprintMature, solver="None", saveat=0.1, callback="None"
)
    ode_prob = ODEProblem(odeFunc, init, tspan, parameters)

    if callback == "None"
        if solver == "None"
            return solve(ode_prob, saveat=saveat)
        else
            return solve(ode_prob, solver, saveat=saveat)
        end
    else
        if solver == "None"
            return solve(ode_prob, saveat=saveat, callback=callback)
        else
            return solve(ode_prob, solver, saveat=saveat, callback=callback)
        end
    end
end

function simulate(
    parameters::ModelParameters, initial_cond::Vector{Float64}, tspan::Tuple{Float64,Float64}, ag_map::AntigenicMap, evorisk_map::EvoRiskMap;
    odeFunc::Function=glOpqua.multistrainImprintMature, #mutateFunc::Function=glOpqua.mutateFunc,
    start_freq::Float64=0.01, solver="None", saveat=0.1
)
    p = copyParams(parameters)
    init = deepcopy(initial_cond)
    ag = copyAntigenicMap(ag_map)
    evorisk = copyEvoRiskMap(evorisk_map)
    t = tspan[1]
    sols = []
    clusters = []
    pars = []

    while t < tspan[2]
        println(("BEGIN", t))
        evorisk_threshold = 1.0 * randexp()
        println(evorisk_threshold)

        # init[end-sum(p.n_strains)+1:end] = -1.0 * evorisk_thresholds
        start_evos = init[end-sum(p.n_strains)+1:end]

        function mutateFunc(u, t, integrator)
            sum(u[end-sum(p.n_strains)+1:end].-start_evos) > evorisk_threshold
            # any(u[end-sum(p.n_strains)+1:end] .> 0.0)
        end

        callback = DiscreteCallback(mutateFunc, terminate!) #TODO: add callback for strain extinction

        push!(sols, simulateODE(p, init, (t,tspan[2]), odeFunc=odeFunc, solver=solver, saveat=saveat, callback=callback))
        push!(clusters, deepcopy(ag.positions))
        push!(pars, copyParams(p))
        t = sols[end].t[end]

        parent_strain, parent_cluster, parent_cluster_num, infecting_host_immunities, infecting_host_immunities_id, same_clusters = sampleInfection(sols[end],p,ag)#TODO:
        # infecting_host_immunities = sampleImmunities(sols[end], parent_strain)#TODO:
        # same_clusters = sameClusters(sols[end])#TODO:

        mutant_pos, cluster_num = mutantFromParent(parent_strain, parent_cluster, parent_cluster_num, infecting_host_immunities, same_clusters, distanceFunc=p.mutantAgDistance) #TODO: make this function wrap mutantAgFromParent, make a function mutantStrainFromParent that makes new strains and also consider that option in wrapper

        if cluster_num < 0 # if new cluster
            ag = glOpqua.antigenicMap( #TODO: Package this into a function that adds a new cluster and creates the model params?
                positions=vcat(ag.positions, mutant_pos'),
                max_cross_immunity=ag.max_cross_immunity,
                K_cross_immunity=ag.K_cross_immunity,
                h_cross_immunity=ag.h_cross_immunity,
                max_imprinted=ag.max_imprinted,
                K_imprinted=ag.K_imprinted,
                h_imprinted=ag.h_imprinted,
                max_matured=ag.max_matured,
                K_matured=ag.K_matured,
                h_matured=ag.h_matured,)
            evorisk = glOpqua.evoRiskMap(
                ag_distances=ag.distances,
                K_evorisk=evorisk.K_evorisk,
                h_evorisk=evorisk.K_evorisk,
                evorisk_by_immunity=hcat(evorisk.evorisk_by_immunity, p.evoRiskMutant(mutant_pos,parent_strain,p)),#TODO: change evorisk_by_immunity to have strains on rows??
                strains_per_cluster=[p.n_strains..., 1])
            p = glOpqua.parameters(
                n_ag_clusters=p.n_ag_clusters + 1,
                n_strains=[p.n_strains..., 1],
                transmission_rates=[p.transmission_rates..., [p.transmissionMutant(mutant_pos,parent_strain,p)]],
                recovery_rates=[p.recovery_rates..., [p.recoveryMutant(mutant_pos,parent_strain,p)]],
                imprinted_loss_rate=p.imprinted_loss_rate,
                matured_loss_rate=p.matured_loss_rate,
                frac_imprinted=ag.frac_imprinted,
                frac_matured=ag.frac_matured,
                immunity_coef=ag.immunity_coef,
                evorisk_coef=evorisk.evorisk_coef,
                evorisk_by_immunity=evorisk.evorisk_by_immunity,
                mutantAgDistance=p.mutantAgDistance,
                birth_rate_uninfected=p.birth_rate_uninfected,
                birth_rates_infected=[p.birth_rates_infected..., [p.birthMutant(mutant_pos,parent_strain,p)]],
                death_rate_uninfected=p.death_rate_uninfected,
                death_rates_infected=[p.death_rates_infected..., [p.deathMutant(mutant_pos,parent_strain,p)]],)
        else
            #TODO: change architecture so strains have IDs, associated clusters, transmission, recovery, birth, death rates, evorisk coefs
            # evorisk = glOpqua.evoRiskMap( #TODO: Package this into a function that adds a new strain and creates the model params?
            #     ag_distances=ag.distances,
            #     K_evorisk=evorisk.K_evorisk,
            #     h_evorisk=evorisk.K_evorisk,
            #     evorisk_by_immunity=hcat(evorisk.evorisk_by_immunity, evoriskMutant(mutant)),#TODO: add evoriskMutant, modify this to work with new strain
            #     strains_per_cluster=p.n_strains + 1)
            # p = glOpqua.parameters(
            #     n_ag_clusters=p.n_ag_clusters,
            #     n_strains=[p.n_strains..., 1], #TODO: modify this to work with new strain
            #     transmission_rates=[p.transmission_rates..., evorisk.transmissionMutant(mutant)],#TODO: add transmissionMutant
            #     recovery_rates=[p.recovery_rates..., evorisk.recoveryMutant(mutant)],#TODO: add recoveryMutant
            #     imprinted_loss_rate=p.imprinted_loss_rate,
            #     matured_loss_rate=p.matured_loss_rate,
            #     frac_imprinted=ag.frac_imprinted,
            #     frac_matured=ag.frac_matured,
            #     immunity_coef=ag.immunity_coef,
            #     evorisk_coef=evorisk.evorisk_coef,
            #     mutantAgDistance=p.mutantAgDistance,
            #     birth_rate_uninfected=p.birth_rate_uninfected,
            #     birth_rates_infected=[p.transmission_rates..., evorisk.birthMutant(mutant)],#TODO: add birthMutant
            #     death_rate_uninfected=p.death_rate_uninfected,
            #     death_rates_infected=[p.death_rates_infected..., evorisk.deathMutant(mutant)],)#TODO: add deathMutant
            println("not yet!")
        end

        init = zeros(p.n_compartments)
        println(length(pars))
        if length(pars) > 0
            println("IN")
            for old_id in pars[end].immunities_ids
                old_id_idx = pars[end].immunities_dict[old_id]
                new_id = [old_id...,0]
                new_id_idx = p.immunities_dict[new_id]
                # init[new_id_idx] = sols[end].u[end][old_id_idx]
                for s in 0:pars[end].n_ag_clusters
                    init[new_id_idx+(s*p.n_immunities)] = sols[end].u[end][old_id_idx+(s*pars[end].n_immunities)]
                end
                # println((infecting_host_immunities_id,old_id))
                if infecting_host_immunities_id == old_id
                    init[new_id_idx+(p.n_ag_clusters*p.n_immunities)] = min(start_freq,init[new_id_idx]) #p.recovery_rates[1][1] / p.transmission_rates[1][1]
                    init[new_id_idx] -= min(start_freq,init[new_id_idx])
                    println("INFECT!")
                end
            end
        end
        init[end-sum(p.n_strains)+1:end-1] = sols[end].u[end][end-sum(pars[end].n_strains)+1:end]
        # println(init)
        #TODO: assign value to newly infected host compartment with mutant

        # println(pars[100000])
        # start_freq = 0.01
        # strain_dist = 9.0 * 500.0 .^ (-1.0 * collect(1.0:sum(p.n_strains)))
        # init = zeros(p.n_compartments)
        # init[1] = 1.0 - (start_freq * sum(p.n_strains))
        # for i in 1:sum(p.n_strains)
        #     init[i*p.n_immunities+1] = start_freq * strain_dist[i]
        #     init[(sum(p.n_strains)+1)*p.n_immunities+i] = start_freq * p.evorisk_coef[(repeat([NAIVE], outer=p.n_ag_clusters), i)]
        # end
    end

    #TODO: process sols into single object?
    return clusters,sols,pars
end
