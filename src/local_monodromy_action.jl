using HomotopyContinuation, LinearAlgebra, Combinatorics, Distances

module LocalMonodromyCode
@enum codes begin
	success
	success_not_branch
	success_action_trivial
	failure_path_crossing_in_monodromy_lift
	failure_in_fiber_path_tracking
end
end

is_failure(S::LocalMonodromyCode.codes) = (S==LocalMonodromyCode.failure_path_crossing_in_monodromy_lift || S==LocalMonodromyCode.failure_in_fiber_path_tracking)
is_success(S::LocalMonodromyCode.codes) = !is_failure(S)

"""
    compute_local_monodromy_action(system, point, vars = variables(system))

Returns detailed information about the local monodromy action of the system at the given point.

"""
function compute_local_monodromy_action(system,point,vars=variables(system))
	# Apply translation to point = 0
	shifted_system = subs(system,vars=>(vars+point))
	# To shift back when outputting results
	shift_back = pt -> pt+point

	# Select random projection to variety dimension
	ambient_dimension = length(vars)
	projection_dimension = ambient_dimension - length(system)	
	subspace_matrix = svd(rand(ComplexF64,projection_dimension,ambient_dimension)).Vt

	# Set up static random plane projection to branch locus dimension
	# and moving base point
	random_orthonormal_basis_for_base_space_matrix = orthogonal_basis_of_random_subspace(projection_dimension,projection_dimension)	
	# Vanishing locus is transpose(A)*(x-base_point) = 0 
	# Need branch locus of dimension at least 3 to require a static
	# slice
	if projection_dimension <= 2 
		vanishing_matrix_for_static_plane = zeros(Complex{Float64},1,projection_dimension)
		substitution_matrix_for_static_plane = random_orthonormal_basis_for_base_space_matrix
	else 
		vanishing_matrix_for_static_plane = conj_transpose(hcat(random_orthonormal_basis_for_base_space_matrix[:,1:(end-3)],random_orthonormal_basis_for_base_space_matrix[:,end]))
		substitution_matrix_for_static_plane = random_orthonormal_basis_for_base_space_matrix[:,(end-2):(end-1)]
	end
	basepoint_path,basepoint_path_variable = path_from_zero(substitution_matrix_for_static_plane)
	basepoint_at_time_equals_one = to_number.(subs(basepoint_path,basepoint_path_variable=>1))


	# In the curve case, complement of line is 0
	if projection_dimension == 1
		matrix_defining_random_line_vanishing = zeros(Complex{Float64},1,projection_dimension)
	else		
		matrix_defining_random_line_vanishing = conj_transpose(random_orthonormal_basis_for_base_space_matrix[:,1:(end-1)])
	end
	# Line is parameterized by t*column+base_point
	column_defining_random_line = random_orthonormal_basis_for_base_space_matrix[:,end]	
	# Find the "localized points" in the fiber for monodromy loops later:

	# Localization filter
	@var variables_for_localization[1:ambient_dimension]
	fiber_system = System([shifted_system;subspace_matrix*vars - basepoint_path],variables=vars,parameters=[basepoint_path_variable])
	fiber_start_points = solutions(solve(fiber_system;target_parameters=[1],show_progress=false,tracker_options=ESPECIALLY_CONSERVATIVE);only_nonsingular=false)
	fiber_end_points = solve(fiber_system,fiber_start_points;start_parameters=[1],target_parameters=[0],tracker_options=ESPECIALLY_CONSERVATIVE,show_progress=false)
	localized_points_in_fiber = [path.start_solution for path in fiber_end_points if sum(abs.(path.solution)) <= DEFAULT_THRESHOLD]

	if length(localized_points_in_fiber) == 0
		return Dict("fiber_points"=>shift_back.(localized_points_in_fiber),
			"orders"=>[1],
			"orbits"=>["trivial"=>[collect(1:length(localized_points_in_fiber))]],
			"system"=>system,
			"variables"=>vars,
			"point"=>point,
			"slice"=>LinearSubspace(subspace_matrix,basepoint_at_time_equals_one-subspace_matrix*point),
			"status"=>LocalMonodromyCode.success_not_branch)
	end
	# Compute a witness superset for the image of the branch locus of the projection

	# (1) Construct the polynomial system which defines the graph of the projection 
	# {(v,x) | F(v) = 0, x - pi(v) = 0}
	variables_for_variety = vars
	@var variables_for_projection[1:projection_dimension]
	graph_system_equations = [shifted_system;subspace_matrix*vars-variables_for_projection]


	# (2) Find a random line in the base space as a parameterization. Use it to substitute into the projection variables for the graph
	@var substitution_variable
	graph_system_equations_with_line_substituted = subs(graph_system_equations,variables_for_projection => (substitution_variable*column_defining_random_line + basepoint_path))

	# (3) Solve the system obtained by adding rank-vanishing conditions from the Jacobian of the system from (2) with respect to the original variety variables
	jacobian_with_respect_to_system_variables = differentiate(graph_system_equations,variables_for_variety)

	nrow,ncol = size(jacobian_with_respect_to_system_variables)
	row_index_combinations = combinations(1:nrow,ambient_dimension)
	column_index_combinations = combinations(1:ncol,ambient_dimension)

	pairs_of_combinations = Base.Iterators.product(row_index_combinations,column_index_combinations)
	# Looking for rank drop less than projection_dimension
	rank_vanishing_equations = [det(jacobian_with_respect_to_system_variables[row_indices,column_indices]) for (row_indices,column_indices) in pairs_of_combinations]
	graph_system_with_rank_vanishing = [graph_system_equations_with_line_substituted;vcat(rank_vanishing_equations...)]	

	# Now we need to randomize down
	system_for_witness_superset = System(graph_system_with_rank_vanishing,variables=[variables_for_variety;substitution_variable],parameters=[basepoint_path_variable])
	initial_system = System(subs(graph_system_with_rank_vanishing,basepoint_path_variable=>1),variables=[variables_for_variety;substitution_variable])

	initial_superset_results = solve(initial_system;show_progress=false,tracker_options=ESPECIALLY_CONSERVATIVE)
	singular_initial_results = singular(initial_superset_results)
	nonsingular_initial_results = nonsingular(initial_superset_results)

	# Attempt deflation if there are singular solutions
	singular_solution_pairs = []
	if length(singular_initial_results)>0
		deflated_system_for_witness_supers,deflated_start_points = deflate(system_for_witness_superset,solutions(singular_initial_results;only_nonsingular=false),[1])
		deflation_deformation_results = solve(deflated_system_for_witness_supers,
			deflated_start_points;
			start_parameters=[1],
			target_parameters=[0],
			tracker_options=ESPECIALLY_CONSERVATIVE,
			show_progress=false)
		# The first ambient_dimension+1 coordinates correspond to
		# the original variables. The rest are auxiliary
		total_dim = ambient_dimension+1
		singular_solution_pairs = [(path.start_solution[1:total_dim],path.solution[1:total_dim]) for path in deflation_deformation_results]
	end

	nonsingular_solution_pairs = []
	if length(nonsingular_initial_results)>0
		nonsingular_deformation_results = solve(system_for_witness_superset,
			solutions(nonsingular_initial_results);
			start_parameters=[1],
			target_parameters=[0],
			tracker_options=ESPECIALLY_CONSERVATIVE,
			show_progress=false)	
		nonsingular_solution_pairs = [(path.start_solution,path.solution) for path in nonsingular_deformation_results]
	end
	# (4) The image (second factor projection) of the points from (3) are the points in the witness superset for the branch locus. Construct a dimension-1 linear system to represent the line
	# using e.g. an SVD computation to find the null space. 

	# Perform localization by deforming superset and retaining start points of paths that end at 0
	all_deformation_solutions = [singular_solution_pairs;nonsingular_solution_pairs]

	filtered_by_localization_superset_solutions = [path_result for path_result in all_deformation_solutions if sum(abs.(path_result[2])) <= DEFAULT_THRESHOLD]	
	# This case can happen if the monodromy action is trivial.
	# In that case, report it
	if length(filtered_by_localization_superset_solutions) == 0
		return Dict("fiber_points"=>shift_back.(localized_points_in_fiber),
			"orders"=>[1],
			"orbits"=>[],
			"system"=>system,
			"variables"=>vars,
			"point"=>point,
			"slice"=>LinearSubspace(subspace_matrix,basepoint_at_time_equals_one-subspace_matrix*point),
			"status"=>LocalMonodromyCode.success_action_trivial)
	end

	all_deformation_start_solutions = (path->path[1]).(all_deformation_solutions)
	real_versions_of_points = [pt[end] for pt in all_deformation_start_solutions]
	real_versions_of_points = [[real(pt),imag(pt)] for pt in real_versions_of_points]

	points_to_encircle_projected_to_line_as_reals = (path->path[1]).(filtered_by_localization_superset_solutions)
	points_to_encircle_projected_to_line_as_reals = (pt->pt[end]).(points_to_encircle_projected_to_line_as_reals)
	points_to_encircle_projected_to_line_as_reals = (pt->[real(pt),imag(pt)]).(points_to_encircle_projected_to_line_as_reals)

	# Now, to count the actual components, we need to 
	# (1) Construct the generator loops
	# (2) Determine the order of each generator (this consequently determines the orbit of all the only-1-generator elements)
	# (3) Compute the orbit over the appropriate product of generators 	
	distances_between_localized_and_other_points = pairwise(Euclidean(),hcat(points_to_encircle_projected_to_line_as_reals...)',hcat(real_versions_of_points...,[0,0])',dims=1)
	loop_radii = [minimum(filter(element->element>0,distances_between_localized_and_other_points[i,:])) for i in 1:length(points_to_encircle_projected_to_line_as_reals)]
	points_and_distances = zip(points_to_encircle_projected_to_line_as_reals,loop_radii)
	

	generator_paths = [box_path_around_point(pt,dist/2,basepoint_at_time_equals_one,column_defining_random_line) for (pt,dist) in points_and_distances]
	fiber_system_general = System([shifted_system;subspace_matrix*vars - variables_for_projection],variables=vars,parameters=variables_for_projection)	

	orbits_and_codes = (generator->determine_generator_orbit(generator,localized_points_in_fiber,fiber_system_general)).(generator_paths)
	orbits = (orbit_result->orbit_result[1]=>orbit_result[2][1]).(zip(generator_paths,orbits_and_codes))
	status_codes = (orbit_result->orbit_result[2]).(orbits_and_codes)
	monodromy_status = LocalMonodromyCode.success
	if any(code->code==LiftingCode.failure,status_codes)
		monodromy_status = LocalMonodromyCode.failure_path_crossing_in_monodromy_lift
	end
	return Dict("fiber_points"=>shift_back.(localized_points_in_fiber),
		"generators"=>generator_paths,
		"orbits"=>orbits,
		"orders"=>(orbit->length(orbit[2])).(orbits),
		"system"=>system,
		"variables"=>vars,
		"point"=>point,
		"slice"=>LinearSubspace(subspace_matrix,basepoint_at_time_equals_one-subspace_matrix*point),
		"status"=>monodromy_status)
end