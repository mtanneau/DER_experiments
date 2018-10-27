"""
    DeferrableLoad

Deferrable Load
"""
mutable struct DeferrableLoad <: Resource
    index::Int  # Unique index of device
    num_timesteps::Int  # Number of time-steps in horizon
    dt::Float64  # Duration of each time-step

    num_cycles::Int  # Number of cycles
    cycle_begin::Vector{Int}  # Beginning of each cycle
    cycle_end::Vector{Int}    # End of each cycle
    cycle_energy_min::Vector{Float64}  # Minimal energy consumption for each cycle
    cycle_energy_max::Vector{Float64}  # Maximal energy consumption for each cycle

    pwr_min::Vector{Float64}  # Minimum power consumption at each time
    pwr_max::Vector{Float64}  # Maximum power consumption at each time
    binary_flag::Bool  # Indicates whether device is hard or soft on-off
end

function DeferrableLoad(;
    index::Integer=0,
    T::Integer=0,
    dt::Float64=1.0,
    num_cycles::Integer=0,
    cycle_begin::Vector{Int}=Int[],
    cycle_end::Vector{Int}=Int[],
    cycle_energy_min::Vector{Float64}=Float64[],
    cycle_energy_max::Vector{Float64}=Float64[],
    pwr_min::Vector{Float64}=Float64[],
    pwr_max::Vector{Float64}=Float64[],
    binary_flag::Bool=true
)
    # Dimension checks
    dt > 0.0 || throw(DomainError("dt must be positive"))
    T == length(pwr_min) || throw(DimensionMismatch(
        "T=$T but pwr_min has size $(length(pwr_min))"
    ))
    T == length(pwr_max) || throw(DimensionMismatch(
        "T=$T but pwr_max has size $(length(pwr_max))"
    ))
    num_cycles == length(cycle_begin) || throw(DimensionMismatch(
        "N=$num_cycles but cycle_begin has size $(length(cycle_begin))"
    ))
    num_cycles == length(cycle_end) || throw(DimensionMismatch(
        "N=$num_cycles but cycle_end has size $(length(cycle_end))"
    ))
    num_cycles == length(cycle_energy_min) || throw(DimensionMismatch(
        "N=$num_cycles but cycle_energy_max has size $(length(cycle_energy_max))"
    ))
    num_cycles == length(cycle_energy_max) || throw(DimensionMismatch(
        "N=$num_cycles but cycle_energy_max has size $(length(cycle_energy_max))"
    ))

    DeferrableLoad(
        index,
        T, dt,
        num_cycles, cycle_begin, cycle_end,
        cycle_energy_min, cycle_energy_max,
        pwr_min, pwr_max,
        binary_flag
    )
end