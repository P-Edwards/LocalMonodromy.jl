using HomotopyContinuation, LinearAlgebra, Combinatorics, Distances, PrettyTables

DEFAULT_THRESHOLD = 1e-6
ESPECIALLY_CONSERVATIVE = TrackerOptions(;automatic_differentiation=4,parameters=:conservative,max_step_size=1e-4,max_initial_step_size=1e-4,max_steps=100_000)


module LiftingCode
@enum codes begin
	success
	failure	
end
end

function determine_generator_orbit(generator,fiber_points,fiber_lift_system)
	previous_acted_on = (pt->Complex.(pt)).(fiber_points)	
	mapped_fiber_points = previous_acted_on
	orbits = []
	while true				
		points_acted_on = follow_system_along_path(fiber_lift_system,generator,previous_acted_on) 	
		min_norms = [minimum([norm(acted_point-fiber_point) for fiber_point in mapped_fiber_points]) for acted_point in points_acted_on]
		maxmin_norm = maximum(min_norms)
		closest_indices = [argmin([norm(acted_point-fiber_point) for fiber_point in mapped_fiber_points]) for acted_point in points_acted_on]
		if maxmin_norm > DEFAULT_THRESHOLD
			return orbits,LiftingCode.failure
		end
		push!(orbits,closest_indices)
		differences = [norm(points_acted_on[i]-mapped_fiber_points[i]) <= DEFAULT_THRESHOLD for i in 1:length(points_acted_on)]	
		if prod(differences)
			break
		end
		previous_acted_on = points_acted_on
	end
	return orbits,LiftingCode.success
end



function local_witness_sets_from_monodromy_action(monodromy_action_output_data)
	if !haskey(monodromy_action_output_data,"fiber_points")
		return Dict("witness_points"=>[[]])
	end
	fiber = monodromy_action_output_data["fiber_points"]
	number_of_fiber_points = length(fiber)	
	total_fiber_orbits = [Set(i) for i in 1:number_of_fiber_points]	
	# The last action is always trivial, so we need to 
	# skip that one 
	for orbit in monodromy_action_output_data["orbits"]
		# In this case we've already identified that 
		# there's only 1 irreducible component and there's
		# nothing else to do
		if length(total_fiber_orbits) == 1
			break
		end			
		for action in orbit[2]
			fiber_transpositions = zip(1:length(fiber),action)
			for (start,stop) in fiber_transpositions						
				if start==stop
					continue
				end
				start_in = [i for i in 1:length(total_fiber_orbits) if start in total_fiber_orbits[i]][1]
				stop_in =  [i for i in 1:length(total_fiber_orbits) if stop in total_fiber_orbits[i]][1]
				if start_in == stop_in
					continue
				end
				merged_orbit = union(total_fiber_orbits[start_in],total_fiber_orbits[stop_in])				
				deleteat!(total_fiber_orbits,sort([start_in,stop_in]))				
				push!(total_fiber_orbits,merged_orbit)
			end
		end					
	end	
	fiber_points_for_witness_sets = [[fiber[i] for i in orbit_set] for orbit_set in total_fiber_orbits]
	system = monodromy_action_output_data["system"]
	vars = monodromy_action_output_data["variables"]
	point = monodromy_action_output_data["point"]
	slice = monodromy_action_output_data["slice"]

	return NumericalLocalIrreducibleDecomposition(
		Dict((length(vars)-length(system))=>[LocalWitnessSet(System(system,variables=vars),slice,fiber_set,point) for fiber_set in fiber_points_for_witness_sets])
			)
end

"""
    NumericalLocalIrreducibleDecomposition{S <: LocalWitnessSet, P <: Int64}

Data structure for holding numerical local irreducible decompositions.

# Fields:
- `Local_Witness_Sets::Dict{P, Vector{S}}`: Dictionary of form (dimension of components)=>Vector of LocalWitnessSets
"""
struct NumericalLocalIrreducibleDecomposition{S<:LocalWitnessSet,P<:Int64}
	Local_Witness_Sets::Dict{P,Vector{S}}	
end


# Close to verbatim from HomotopyContinuation.jl - numerical_irreducible_decomposition
function max_dim(N::NumericalLocalIrreducibleDecomposition)
	k = keys(N.Local_Witness_Sets)
	if !isempty(k)
        maximum(k)
    else
        -1
    end
end

function Base.show(io::IO, N::NumericalLocalIrreducibleDecomposition)
    D = N.Local_Witness_Sets
    if !isempty(D)
        total = sum(length(last(Ws)) for Ws in D)
    else
        total = 0
    end
    header = "\n Numerical local irreducible decomposition with $total components" 
    println(io, header)

    mdim = max_dim(N)
    if mdim >= 0
        for d = max_dim(N):-1:0
            if haskey(D, d)
                ℓ = length(D[d])
                if ℓ > 0
                    println(io, "• $ℓ component(s) of dimension $d.")
                end
            end
        end

        println(io, "\n degree table of components:")
        degree_table(io, N)
    end
end
function degree_table(io, N::NumericalLocalIrreducibleDecomposition)
    D = N.Local_Witness_Sets
    k = collect(keys(D))
    sort!(k, rev = true)
    n = length(k)

    headers = ["dimension", "degrees of components"]
    data = Matrix{Union{Int,String}}(undef, n, 2)

    for (i, key) in enumerate(k)
        data[i, 1] = key
        components = Tuple(length(W.R) for W in D[key])
        if length(components) == 1
            data[i, 2] = first(components)
        elseif length(components) <= 10
            data[i, 2] = string(components)
        else
            s = string(components[1:10])
            data[i, 2] = string("(", s[2:end-1], ", ...)")
        end
    end

    PrettyTables.pretty_table(
        io,
        data;
        header = headers,
        tf = PrettyTables.tf_unicode_rounded,
        alignment = :c,
        header_crayon = PrettyTables.Crayon(bold = false),
        border_crayon = PrettyTables.Crayon(faint = true),
    )
end

## End verbatim block

"""
    NumericalLocalIrreducibleDecomposition(system, point, vars = variables(system); max_tries = 3)

Compute decomposition of `system` at `point`. Result is of type NumericalLocalIrreducibleDecomposition.

# Example
```julia-repl
julia> using HomotopyContinuation
julia> @var x y z
julia> result = NumericalLocalIrreducibleDecomposition([x^2+y^2+z^2],[0;0;0])
 Numerical local irreducible decomposition with 1 components
• 1 component(s) of dimension 2.

 degree table of components:
╭───────────┬───────────────────────╮
│ dimension │ degrees of components │
├───────────┼───────────────────────┤
│     2     │           2           │
╰───────────┴───────────────────────╯

```

# Arguments:
- `system`: Vector of ModelKit expressions.
- `point`: Vector of complex numbers, 1 for each variable. 
- `vars`: Optional vector for variable order.
- `max_tries`: Number of retries when numerical failure detected.
"""
function NumericalLocalIrreducibleDecomposition(system,point,vars=variables(system);max_tries=3)
	i = 1
	action_data = compute_local_monodromy_action(system,point,vars)
	while i < max_tries && is_failure(action_data["status"])
		action_data = compute_local_monodromy_action(system,point,vars)
		i+=1
	end
	if is_failure(action_data["status"])
		error("Decomposition failed after maximum number of attempts.")
	end
	
	return local_witness_sets_from_monodromy_action(action_data)
end

"""
    nlid(system, point, vars = variables(system))

Alternate name for `NumericalLocalIrreducibleDecomposition`.

"""
function nlid(system,point,vars=variables(system);max_tries=3)
	return NumericalLocalIrreducibleDecomposition(system,point,vars;max_tries=max_tries)
end

function count_local_components(nlid_for_system)
	local_witness_sets = nlid_for_system.Local_Witness_Sets
	dimensions = keys(local_witness_sets)
	top_dim_sets = local_witness_sets[dimensions[1]]
	return [length(loc_witness_set.R) for loc_witness_set in top_dim_sets]
end

"""
    decompose_branch_locus(system, vars = variables(system), codimension = length(system))

Returns a numerical irreducible decomposition of the branch locus of the system.

# Example
```julia-repl
julia> using HomotopyContinuation
julia> @var x y z
julia> result = decompose_branch_locus([x^2+y^2+z^2])
 Numerical irreducible decomposition with 2 components
• 2 component(s) of dimension 1.

 degree table of components:
╭───────────┬───────────────────────╮
│ dimension │ degrees of components │
├───────────┼───────────────────────┤
│     1     │        (1, 1)         │
╰───────────┴───────────────────────╯

```

# Arguments:
- `system`: Vector{ModelKit.Expression}
- `vars`: Vector{ModelKit.Variable}
- `codimension`: Integer
"""
function decompose_branch_locus(system,vars=variables(system),codimension=length(system))
	# Select random projection to variety dimension
	ambient_dimension = length(vars)
	projection_dimension = ambient_dimension - length(system)	
	subspace_matrix = rand(Complex{Float64},projection_dimension,ambient_dimension)

	graph_system_equations = [system;subspace_matrix*vars]
	jacobian_of_system = differentiate(graph_system_equations,vars)
	nrow,ncol = size(jacobian_of_system)
	row_index_combinations = combinations(1:nrow,ambient_dimension)
	column_index_combinations = combinations(1:ncol,ambient_dimension)

	pairs_of_combinations = Base.Iterators.product(row_index_combinations,column_index_combinations)
	# Looking for rank drop less than projection_dimension
	rank_vanishing_equations = [det(jacobian_of_system[row_indices,column_indices]) for (row_indices,column_indices) in pairs_of_combinations]
	graph_system_with_rank_vanishing = [system;vcat(rank_vanishing_equations...)]	


	system_with_rank = System(graph_system_with_rank_vanishing,variables=vars)
	# The numerical irreducible decomposition function doesn't play nice with
	# non-complete intersections at the moment. So we'll randomize down manually
	# and filter out the witness sets at the end which don't correspond to 
	# components of the original system	
	with_extra = nid(system_with_rank)
	for (dim,this_dim_sets) in with_extra.Witness_Sets 		
		to_delete = []
		for index in length(this_dim_sets)
			ws = this_dim_sets[index]
			first_solution = solutions(ws)[1]
			# With probability 1 the extra components introduced by randomization 
			# don't intersect the originals, and so we can expect to see only 
			# witness points that don't satisfy the original equations
			how_close_to_satisfying_original_system = norm(to_number(sum(subs(system_with_rank.expressions,vars=>first_solution))))
			if how_close_to_satisfying_original_system > DEFAULT_THRESHOLD
				push!(to_delete,index)
			end
		end
		if length(to_delete)==length(this_dim_sets)
			delete!(with_extra.Witness_Sets,dim)
		else
			deleteat!(this_dim_sets,to_delete)	
		end	
	end
	return with_extra
end



