# Fixed load

## Parameters

* $T$: number of time-steps
* $(P_{1}, ..., P_{T})$: fixed load at each time step

## Variables
* $(p_{1}, ..., p_{T})$: net load at each time step

## Model

$$
    \begin{array}{ccl}
        \displaystyle \min_{p} & 0\\
        s.t. & p_{t} &= P_{t} \ \ \forall t
    \end{array}
$$


# Curtailable load

## Parameters

* $T$: numver of time-steps
* $(P_{1}, ..., P_{T})$: power consumption at time $t$, if device is on

## Variables

* $p_{t}$: net load at time $t$
* $u_{t}$: binary indicator if device is on at time $t$

## Model

$$
    \begin{array}{ccl}
        \displaystyle \min_{p, u} & 0\\
        s.t. 
        & p_{t} &= u_{t} P_{t} \ \ \forall t\\
        & u_{t} & \in \mathbb{B} \ \ \forall t
    \end{array}
$$


# Uninterruptible load

## Parameters

* $T$: numver of time-steps
* $L$: cycle's duration
* $(P_{1}, ..., P_{L})$: power consumption at step $l$ of the cycle

## Variables

* $p_{t}$: net load at time $t$
* $v_{t}$: binary, indicates whether cycle was started at time $t$

## Model

$$
    \begin{array}{ccl}
        \displaystyle \min_{p, u} & 0\\
        s.t. 
        & p_{t} - \displaystyle \sum_{l=1}^{L} P_{l} v_{t+1-l}  &= 0 \ \forall t\\
        & \displaystyle \sum_{t=1}^{T} v_{t} &=1\\
        & v_{t} & \in \mathbb{B} \ \ \forall t
    \end{array}
$$

# Deferrable load

## Parameters
* $P^{min}_{t}$: minimum power at time $t$
* $P^{max}_{t}$: maximum power at time $t$
* $E^{min}$: minimum energy requirement over horizon
* $E^{max}$: maximum energy requirement over horizon

## Variables
* $p_{t}$: net load at time $t$
* $u_{t}$: on-off indicator at time $t$

## Model

$$
    \begin{array}{ccl}
        \displaystyle \min_{p, u} & 0\\
        s.t. 
        & E^{min} \leq \sum_{t} p_{t} &\leq E^{max}\\
        & u_{t}  P^{min}_{t} \leq p_{t} &\leq u_{t} P^{max}_{t} \ \ \forall t\\
        & u_{t} & \in \mathbb{B} \ \ \forall t
    \end{array}
$$


# Battery

## Parameters
* $P^{chg, min}_{t}$: Minimum charging power at time $t$
* $P^{chg, max}_{t}$: Maximum charging power at time $t$
* $P^{dis, min}_{t}$: Minimum discharging power at time $t$
* $P^{dis, max}_{t}$: Maximum discharging power at time $t$
* $E^{min, max}$: Minimum and maximum state-of-charge of the battery
* $η^{chg, dis}$: Charging/Discharging efficiencies of the battery

## Variables
* $p_{t}$: Net load at time $t$
* $p^{chg}_{t}, p^{dis}_{t}$: Charging/Discharging power at time $t$
* $e_{t}$: State-of-charge at time $t$
* $u^{chg}_{t}, u^{dis}_{t}$: Charging/Discharging indicator at time $t$

## Model

$$
    \begin{array}{ccl}
        \displaystyle \min_{p, u, e} & 0\\
        s.t.
        & p_{t} &= p^{chg}_{t} - p^{dis}_{t} \ \ \forall t\\ 
        & E^{min} \leq e_{t} &\leq E^{max}  \ \ \forall t\\
        & e_{t} - e_{t-1} &= \eta^{chg} p^{chg}_{t} - \dfrac{1}{\eta^{dis}} p^{dis}_{t} \ \ \forall t\\
        & u^{chg}_{t}  P^{chg, min}_{t} \leq p^{chg}_{t} &\leq u^{chg}_{t} P^{chg, max}_{t} \ \ \forall t \\
        & u^{dis}_{t}  P^{dis, min}_{t} \leq p^{dis}_{t} &\leq u^{dis}_{t} P^{dis, max}_{t} \ \ \forall t \\
        & u^{chg}_{t} + u^{dis}_{t} & \leq 1 \ \ \forall t\\
        & u^{chg}_{t}, u^{dis}_{t} & \in \mathbb{B} \ \ \forall t
    \end{array}
$$

# Thermal load

## Parameters
* $P^{min}_{t}, P^{max}_{t}$: Minimum/maximum power at time $t$
* $\Theta^{min}_{t}, \Theta^{max}_{t}$: Minimum/maximum inside temperature at time $t$
* $\Theta^{ext}_{t}$: Outside temperature at time $t$
* $η$: Thermal efficiency
* $\mu$: Thermal conductance
* $C$: Thermal capacity

## Variables
* $p_{t}$: Power consumption at time $t$
* $\theta_{t}$: Inside temperature at time $t$
* $u_{t}$: On-off indicator at time $t$

## Model

$$
    \begin{array}{cl}
        \displaystyle \min_{θ, p, u} & 0\\
        s.t.
        & Θ^{min}_{t} \leq θ_{t} \leq Θ^{max}_{t} &\forall t\\
        & \dfrac{\mu}{C} (\Theta^{ext}_{t} - \theta_{t}) + \dfrac{\eta}{C} p_{t} = \dfrac{1}{\Delta \tau}(\theta_{t+1} - \theta_{t}) &\forall t\\
        & u_{t}  P^{min}_{t} \leq p_{t} \leq u_{t} P^{max}_{t} &\forall t \\
        & u_{t}  \in \mathbb{B} &\forall t
    \end{array}
$$