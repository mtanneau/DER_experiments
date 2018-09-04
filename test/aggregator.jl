import MathProgBase
const MPB = MathProgBase

T = 2
netload_min = 0.0 * ones(T)
netload_max = 10.0 * ones(T)
p_mkt = [2.0, 3.0]
p_tou = [1.0, 3.0]

# build houses
b = DR.DER.Battery(
    T=T,
    dt=1.0,
    soc_min=zeros(T), soc_max=10.0*ones(T), soc_init=0.0, selfdischargerate=1.0,
    charge_pwr_min=0.0, charge_pwr_max=6.0, discharge_pwr_min=0.0, discharge_pwr_max=6.0,
    charge_eff=1.0, discharge_eff=1.0
)
l1 = DR.DER.FixedLoad(T=T, load=ones(T))
l2 = DR.DER.FixedLoad(T=T, load=2.0*ones(T))
h1 = DR.DER.House(T=T, netload_min=zeros(T),netload_max=5.0*ones(T), price=0.0*p_tou, appliances=[b, l1])
h2 = DR.DER.House(T=T, netload_min=zeros(T),netload_max=5.0*ones(T), price=0.0*p_tou, appliances=[b, l2])

# Build aggregator
a = DR.Aggregator(T, netload_min, netload_max, p_mkt, [h1, h2])

# Build centralized model
m = DR.build_centralized_model(a, CplexSolver(CPX_PARAM_SCRIND=0))
MPB.optimize!(m)

@printf("Obj: %f\n", MPB.getobjval(m))

# Solve DW relaxation

