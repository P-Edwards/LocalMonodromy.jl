module LocalMonodromy

export NumericalLocalIrreducibleDecomposition
export nlid
export compute_local_monodromy_action
export decompose_branch_locus


include("utilities.jl")
include("local_witness_set.jl")
include("local_monodromy_action.jl")
include("nlid.jl")


end