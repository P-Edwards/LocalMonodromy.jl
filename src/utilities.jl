using HomotopyContinuation, LinearAlgebra

# Return a matrix populated by a random orthonormal basis, complex
function orthogonal_basis_of_random_subspace(ambient_dimension,projection_dimension)
	@assert projection_dimension <= ambient_dimension
	random_matrix = rand(Complex{Float64},ambient_dimension,projection_dimension)
	left_singular_vectors = svd(random_matrix).U
	return left_singular_vectors[:,1:projection_dimension]
end

function path_from_zero(matrix_for_static_plane)
	@unique_var t
	random_point = rand(Complex{Float64},size(matrix_for_static_plane,2))
	MAX_RADIUS = 0.1
	max_norm_of_rand = norm([Complex(1.0,1.0) for _ in length(random_point)])
	random_point = (MAX_RADIUS/max_norm_of_rand)*random_point
	random_point = matrix_for_static_plane*random_point	
	return (t*random_point,t)
end

function conj_transpose(input_matrix)
	return conj(transpose(input_matrix))
end


function box_path_around_point(point_to_encircle,max_radius,point_to_start_from,column_in_ambient)
	point_to_start_from = Complex.(point_to_start_from)
	path_in_real = [point_to_encircle + max_radius*[0;-1],point_to_encircle+max_radius*[-1;0],point_to_encircle+max_radius*[0;1],point_to_encircle+max_radius*[1;0],point_to_encircle+max_radius*[0;-1]]
	path_in_complex_line = (pt->Complex(pt[1],pt[2])).(path_in_real)
	path_in_ambient_space = (pt->pt*column_in_ambient+point_to_start_from).(path_in_complex_line)
	return [[point_to_start_from];path_in_ambient_space;[point_to_start_from]]
end

function follow_system_along_path(input_system,path_to_follow,points_to_follow)
	current_points = points_to_follow
	path_edges = [(path_to_follow[i],path_to_follow[i+1]) for i in 1:(length(path_to_follow)-1)]
	max_condition = 0
	for edge in path_edges
		straight_line_homotopy = ParameterHomotopy(input_system;start_parameters=edge[1],target_parameters=edge[2])
		homotopy_results = solve(straight_line_homotopy,current_points;show_progress=false,tracker_options=ESPECIALLY_CONSERVATIVE)
		current_points = solutions(homotopy_results)			
	end
	return current_points
end

function rank_of_expression_matrix(mat;expect_deficiency=false)
	# Coercing ModelKit to output an actual numerical 
	# matrix that LinearAlgebra will accept takes a 
	# bit of a cludge. Otherwise, the entries can
	# be of type ModelKit.Expression.
	calculated_rank = rank(Complex.(to_number.(expand.(mat))))
	# The  numerical rank computation here occasionally disagrees with
	# the rank computed during tracking. At the risk of degrading 
	# the generality of this function, the following handles
	# a common case directly. 
	if expect_deficiency && calculated_rank == minimum(size(mat))
		return minimum(size(mat))-1
	else
		return calculated_rank
	end
end


function deflate(system,points,initial_parameters;current_tries=0,max_tries=3)
	# Stopping criterion if deflation doesn't seem 
	# to be working
	if current_tries == max_tries
		return system,points
	end

	# Basic values	
	params = system.parameters
	vars = system.variables
	number_of_variables = length(vars)
	number_of_functions = length(system.expressions)
	system_at_initial = subs(system.expressions,params=>initial_parameters)
	jacobian_at_initial = differentiate(system_at_initial,vars)
	rank_at_points = rank_of_expression_matrix(subs(jacobian_at_initial,vars=>points[1]),expect_deficiency=true)
	if rank_at_points == number_of_functions
		rank_at_points = number_of_functions-1
	end
	dimnull = number_of_variables - rank_at_points

	# Set up deflation system
	random_vectors = rand(ComplexF64,number_of_variables,number_of_variables)
	@unique_var lambdas[1:rank_at_points,1:dimnull]
	identity_matrix = Matrix(Complex(1,0)I,dimnull,dimnull)	
	inflated_variables = [vars;vcat(lambdas...)]
	new_parameterized_system = System(vcat(system.expressions...,differentiate(system.expressions,vars)*random_vectors*vcat(identity_matrix,lambdas)...),variables=inflated_variables,parameters=params)

	# Compute new start solutions for this system.
	new_system_at_initial_expressions = subs(new_parameterized_system.expressions,params=>initial_parameters)
	new_points = []
	enough_rank = false
	for point in points 
		matrix_at_this_point = Complex.(to_number.(expand.(subs(jacobian_at_initial,vars=>point)*random_vectors)))
		qone = matrix_at_this_point[:,1:dimnull]
		qtwo = matrix_at_this_point[:,(dimnull+1):end]
		solutions_matrix = Matrix{Any}(undef,number_of_variables-dimnull,dimnull)
		# This uses Julia's base LinearAlgebra, which can lose precision. 
		for index in 1:dimnull
			solutions_matrix[:,index] =  qtwo \ -qone[:,index]
		end
		point_with_multipliers = [point;vcat(solutions_matrix...)]
		push!(new_points,point_with_multipliers)
		jacobian_of_deflated_system_at_this_point = subs(differentiate(new_system_at_initial_expressions,inflated_variables),inflated_variables=>point_with_multipliers)
		rank_deflated_enough_at_this_point = (rank_of_expression_matrix(jacobian_of_deflated_system_at_this_point) == length(inflated_variables))
		enough_rank = enough_rank || rank_deflated_enough_at_this_point
	end
	# Several deflation steps may be necessary
	if enough_rank				
		return new_parameterized_system,new_points
	else
		# In this case we can square up the system
		number_new_functions = length(new_parameterized_system.expressions)
		random_down_matrix = rand(ComplexF64,length(variables(new_parameterized_system)),number_new_functions)
		squared_up_parameterized = System(random_down_matrix*new_parameterized_system.expressions,variables=variables(new_parameterized_system),parameters=params)
		return deflate(squared_up_parameterized,new_points,initial_parameters;current_tries=current_tries+1,max_tries=max_tries)
	end
end