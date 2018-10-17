T = 2
netload_min = 0.0 * ones(T)
netload_max = 10.0 * ones(T)
p_mkt = [2.0, 3.0]
p_tou = [1.0, 3.0]

# build houses
b1 = DR.DER.Battery(
    index=1,
    T=T,
    dt=1.0,
    soc_min=zeros(T), soc_max=10.0*ones(T), soc_init=0.0, selfdischargerate=1.0,
    charge_pwr_min=0.0, charge_pwr_max=6.0, discharge_pwr_min=0.0, discharge_pwr_max=6.0,
    charge_eff=1.0, discharge_eff=1.0
)
b2 = DR.DER.Battery(
    index=2,
    T=T,
    dt=1.0,
    soc_min=zeros(T), soc_max=10.0*ones(T), soc_init=0.0, selfdischargerate=1.0,
    charge_pwr_min=0.0, charge_pwr_max=6.0, discharge_pwr_min=0.0, discharge_pwr_max=6.0,
    charge_eff=1.0, discharge_eff=1.0
)
l1 = DR.DER.FixedLoad(index=1, T=T, load=ones(T))
l2 = DR.DER.FixedLoad(index=2, T=T, load=2.0*ones(T))
h1 = DR.DER.House(index=1, T=T, netload_min=zeros(T),netload_max=5.0*ones(T), price=0.0*p_tou, appliances=[b1, l1])
h2 = DR.DER.House(index=2, T=T, netload_min=zeros(T),netload_max=5.0*ones(T), price=0.0*p_tou, appliances=[b2, l2])

# Build aggregator
agg = DR.Aggregator(T, netload_min, netload_max, p_mkt, [h1, h2])

# Build centralized model
m = DR.buildmodel(agg, GLPKSolverMIP())
MPB.optimize!(m)

@printf("Obj: %f\n", MPB.getobjval(m))

# Solve DW relaxation
