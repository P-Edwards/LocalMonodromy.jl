using HomotopyContinuation
export LocalWitnessSet


AbstractSystem = HomotopyContinuation.AbstractSystem
AbstractSubspace = HomotopyContinuation.AbstractSubspace

"""
    LocalWitnessSet

Data type for storing local witness sets

# Fields:
- `F::AbstractSystem`: System
- `L::AbstractSubspace`: Linear slice
- `R::Vector{Vector{Complex64}}`: Localized witness points
- `projective::Bool`: Always false. For compatibility with HomotopyContinuation.jl
- `point::Vector{Complex64}`: Point at which localized.
"""
struct LocalWitnessSet{
	S<:AbstractSystem,
	Sub<:AbstractSubspace,
	R<:Union{Vector{Complex},PathResult},
	P<:Vector{Complex}
	}
	F::S
    L::Sub    
    R::Vector{R}
    projective::Bool
    point::P
end

"""
    LocalWitnessSet(F::Union{Vector{Expression}, System}, L::LinearSubspace, R::Union{Vector{Number}, Vector{Complex}}, point::Union{Vector{Number}, Vector{Complex}})

Initialize a LocalWitnessSet. 

# Arguments:
- `F`: Polynomial system.
- `L`: Linear subspace.
- `R`: Localized witness points.
- `point`: Point at which localized.
"""
function LocalWitnessSet(
    F,
    L::LinearSubspace,
    R,
    point    
)
	f = fixed(System(F),compile=false)
	r = (vec->Vector{Complex}(vec)).(R)
	pt = Vector{Complex}(point)
    LocalWitnessSet(f, L, r,false,pt)
end