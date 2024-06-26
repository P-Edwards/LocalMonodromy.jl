LocalMonodromy.jl
=======================================

Version
-------
1.0.0


Copyright (C) 2024 [Parker
Edwards](https://parkeredw.com)


Installation
-------------
For global installation from Julia in package mode (press `]` to enter package mode, backspace to exit): `add https://github.com/P-Edwards/LocalMonodromy.jl.git`

To install as a standalone local package: Clone this repository to a directory. Then, from Julia in package mode: `activate <path/to/cloned/copy/of/project/root>` followed by `instantiate`. 


Basic usage
------
To use multiple cores, run `export JULIA_NUM_THREADS=<number_to_use>` before starting Julia. Examples are available at [https://github.com/P-Edwards/local-monodromy-examples](https://github.com/P-Edwards/local-monodromy-examples).


This package primarily exposes the following function: 

	NumericalLocalIrreducibleDecomposition(system,point,vars=variables(system);max_tries=3)

Where 
* `system` is a vector of polynomial expressions
* `point` is a `Vector{Complex}` at which to compute the numerical local irreducible decomposition.
* `vars` is a `Vector{ModelKit.Variable}` specifying the variable order.
* `max_tries` is an `Int` specifying how many times to retry if numerical failures are detected.

The output is a `NumericalLocalIrreducibleDecomposition`, a type which is overviewed below.

**Example**: Compute a numerical local irreducible decomposition.
	
	julia> using LocalMonodromy, HomotopyContinuation
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





The output `result` is of type `NumericalLocalIrreducibleDecomposition` which has one field:
* `result.Local_Witness_Sets` is a `Dict(Int,Vector{LocalWitnessSets})`, which contains 1 vector of `LocalWitnessSet` for each dimension. Currently, only pure dimension is supported. 

**Example**

	julia> result.Local_Witness_Sets[2]
	1-element Vector{LocalWitnessSet}

Each `LocalWitnessSet` has the following fields: 
* `F` the polynomial system
* `L` a `LinearSubspace` for slicing
* `R` a vector of localized witness points in `V(F)\cap R` 
* `point` the point at which the input was localized



License
-------
This repository is distributed under GPLv3. 
