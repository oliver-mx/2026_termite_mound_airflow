# 24h termite mound air flow simulation with optimal radius and height.

using Pkg; Pkg.activate("."); Pkg.instantiate();
using Trixi, OrdinaryDiffEqLowStorageRK, Interpolations, QuadGK, FastGaussQuadrature, Plots
using Trixi: AbstractEquations, @muladd
import Interpolations: Line
import Trixi: flux_ranocha, ln_mean, inv_ln_mean, flux, varnames, cons2cons, cons2prim, prim2cons, cons2entropy, max_abs_speeds

include("src/equations/termite_mound_1d.jl")
include("src/callbacks/update_velocity_plot_24h.jl")

###############################################################################
# Output

out_dir = joinpath(@__DIR__, "out")
if isdir(out_dir) == false
    mkdir(out_dir)
end
output_dir = joinpath(out_dir, "reference_mound")
if isdir(output_dir) == false
    mkdir(output_dir)
end

@inline function f_main(r,h)
###############################################################################
# Semidiscretization

(γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr², ρₕ₀, T_ref, t_ref, T0, v0, Ti_LI) = termite_parameters(r,h)
equations = TermiteMoundEquations1D(;γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr², ρₕ₀, T_ref, t_ref, T0, v0, Ti_LI)

initial_condition = TermiteMoundInitialCondition
volume_flux = flux_ranocha
surface_flux = flux_ranocha

dg = DGSEM(polydeg = 4, surface_flux = flux_ranocha, volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

mesh = TreeMesh((0.0,), (1.0,), initial_refinement_level = 5, n_cells_max = 1000, periodicity = true) 

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg, source_terms = source_terms, boundary_conditions = boundary_condition_periodic);

###############################################################################
# ODE solvers, callbacks etc.

tspan = (-0.1, 1.0).* 86400 ./ equations.tᵣ

ode = semidiscretize(semi, tspan)

stepsize_callback = StepsizeCallback(cfl = 0.4)
update_velocity_callback = UpdateVelocityCallback(CarpenterKennedy2N54(williamson_condition = false), output_dir)

amr_controller = ControllerThreeLevel(semi, IndicatorMax(semi, variable = (u, equations) -> u[2]), base_level = 5, max_level = 6, max_threshold = -0.35)
amr_callback = AMRCallback(semi, amr_controller, interval = 1, adapt_initial_condition = true, adapt_initial_condition_only_refine = true)

callbacks = CallbackSet(stepsize_callback, update_velocity_callback, amr_callback)

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false); dt = 1.0,
            ode_default_options()..., callback = callbacks);

    return (callbacks.discrete_callbacks[2].affect!.e(), callbacks.discrete_callbacks[2].affect!.b())
end

###############################################################################
# Simulation

t_start = time()
(ϵ, ζ) = f_main(0.6,2.0)
t_end = time() - t_start

println("(ϵ, ζ) = ($(round(ϵ, digits=4)), $ζ)")
println("Simulation time: $(round(t_end,digits=1))s")

include("src/plots/reference_mound_plot.jl")
