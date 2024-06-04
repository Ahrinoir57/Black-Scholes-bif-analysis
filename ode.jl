using Revise
using Test, ForwardDiff, Parameters
using Plots
# using GLMakie; Makie.inline!(true)
using BifurcationKit, Test
const BK = BifurcationKit
####################################################################################################

# vector field
function BlackScholes(z, p, t = 0)
	@unpack r, s, σ, c = p
	S, dS = z
	[
        dS,
	    2.0 * (c * S - r * s * dS) / (σ^2 * s^2)
	]
end


# parameter values
di = Dict{String, Float64}()

f = open("julia_scripts/julia_data.txt")
z0 = [0., 0.]

for line in eachline(f)
	if line == "----"
		
		par_tm = (r = di["r"], c = di["c"], s = di["s"], σ = di["σ"])
		# Bifurcation Problem
		prob = BifurcationProblem(BlackScholes, z0, par_tm, (@lens _.r);
			record_from_solution = (x, p) -> (S = x[1], dS = x[2]),)

		# continuation parameters
		bif_par = "r"
		optcont = ContinuationPar(nev= 5, max_steps = 100,
			p_min = 0.5 * di[bif_par], p_max = 2 * di[bif_par], ds = 0.01 * di[bif_par], dsmin =  0.001 * di[bif_par], dsmax = di[bif_par], tol_stability = 2e-6)

		# compute the branch of solutions
		br = continuation(prob, PALC(tangent=Bordered()), optcont; bothside = true, normC = norminf)

		open("./julia_scripts/output.txt", "a") do file
			write(file, "Params:\n")
			for (k,v) in zip(keys(par_tm), par_tm)
				write(file, repr(k) * "#" * repr(v) * "\n")
			end
			write(file, "Result:\n")
			write(file, repr(br) * "\n")
			for i in 1:14
				write(file, repr(br[i]) * "\n")
			end
			write(file, "-------\n")
		end

	else
		var, val = split(line, "=")
		if strip(var) == "sigma"
			di["σ"] = parse(Float64, val)
		elseif strip(var) == "S"
			z0[1] = parse(Float64, val)
		elseif strip(var) == "dS"
			z0[2] = parse(Float64, val)
		else
			di[strip(var)] = parse(Float64, val)
		end
	end
end
