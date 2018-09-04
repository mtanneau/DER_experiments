# Build an operational model of a house with each appliance
import Cbc:CbcSolver
import CPLEX:CplexSolver
import Linda

T = 2  # two time-periods
p = [0.0, 1.0]  # Price vector

appliances = []

b = DR.DER.Battery(
    T=T,
    soc_min=zeros(T), soc_max=ones(T),
    charge_pwr_max=1.0, discharge_pwr_max=1.0,
    charge_eff=1.0, discharge_eff=1.0
)
push!(appliances, b)

l = DR.DER.FixedLoad(
    T=T, load=ones(T)
)
push!(appliances, l)

h = DR.DER.House(0, T, 1.0, -5*ones(T), 5*ones(T), p, [b, l])

env = Linda.LindaEnv()
o = Linda.Oracle.LindaOracleMIP(h, CplexSolver(CPX_PARAM_SCRIND=0))
Linda.Oracle.query!(o, zeros(2*T), 0.0)

@test Linda.Oracle.get_sp_dual_bound(o) ≈ 0.0