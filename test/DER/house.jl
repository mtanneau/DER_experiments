# Build an operational model of a house with each appliance

T = 2  # two time-periods
p = [0.0, 1.0]  # Price vector

appliances = []

battery = DR.DER.Battery(
    T=T,
    soc_min=zeros(T), soc_max=ones(T),
    charge_pwr_max=1.0, discharge_pwr_max=1.0,
    charge_eff=1.0, discharge_eff=1.0
)
push!(appliances, battery)

fixedLoad = DR.DER.FixedLoad(T=T, load=ones(T))
push!(appliances, fixedLoad)

curtLoad = DR.DER.CurtailableLoad(T=T, load=ones(T))
push!(appliances, curtLoad)

deferLoad = DR.DER.DeferrableLoad(
    T=T,
    num_cycles=1, cycle_begin=[1], cycle_end=[T],
    cycle_energy_min=[1.0], cycle_energy_max=[1.0],
    pwr_min=zeros(T), pwr_max=ones(T),
    binary_flag=false
)
push!(appliances, deferLoad)

h = DR.DER.House(0, T, 1.0, -5*ones(T), 5*ones(T), p, appliances)

env = Linda.LindaEnv()
o = Linda.Oracle.LindaOracleMIP(h, GLPKSolverMIP())
Linda.Oracle.query!(o, zeros(2*T), 0.0)

@test length(Linda.Oracle.get_new_columns(o)) >= 1