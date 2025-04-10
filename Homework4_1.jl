# Define constants
using Plots
T = 700 + 273  # Temperature in Kelvin
R = 8.314  # Universal gas constant, J/(mol·K)
F = 96485  # Faraday's constant, C/mol
n = 8  # Number of electrons transferred

# Partial pressures (atm)
P_O2 = 0.18  # Oxygen in cathode
P_CH4 = 0.60  # Methane in anode
P_CO2 = 0.20  # Carbon dioxide in anode
P_H2O = 0.20  # Water vapor in anode

# Standard cell potential (V)
E_0 = 1.1

# Compute Nernst potential
E_Nernst = E_0 + (R * T / (n * F)) * log((P_CO2 * P_H2O^2) / (P_CH4 * P_O2^2))

# Exchange current density equations
function i0_cathode(T, C_O2)
    return 3.8e6 * exp(-8170 / T) * C_O2
end

function i0_anode(T, C_CH4)
    return 1.3e7 * exp(-8427 / T) * C_CH4
end

# Compute concentrations using ideal gas law (C = P / (RT))
C_O2 = P_O2 * 101325 / (R * T)
C_CH4 = P_CH4 * 101325 / (R * T)

# Compute exchange current densities
i0_c = i0_cathode(T, C_O2)
i0_a = i0_anode(T, C_CH4)

# Given current densities in A/cm² (convert to A/m²)
j_values = [i * 1e4 for i in 0.1:0.1:3.0]  # Range of current densities for plotting

# Charge transfer coefficient
alpha = 0.5

# Activation overpotential function
function activation_overpotential(j, i0)
    return (R * T / (alpha * n * F)) * asinh(j / (2 * i0))
end

# Compute activation overpotentials
eta_act_anode = [activation_overpotential(j, i0_a) for j in j_values]
eta_act_cathode = [activation_overpotential(j, i0_c) for j in j_values]

# Diffusion properties
d_O2 = 3.66e-7  # m²/s
L_cathode = 10e-6  # m
d_CH4 = 9.66e-7  # m²/s
L_anode = 200e-6  # m

# Limiting current density calculation
j_lim_cathode = (n * F * d_O2 * C_O2) / L_cathode
j_lim_anode = (n * F * d_CH4 * C_CH4) / L_anode

# Concentration overpotential function
function concentration_overpotential(j, j_lim)
    return (R * T / (n * F)) * log(j_lim / (j_lim - j))
end

# Compute concentration overpotentials
eta_conc_anode = [concentration_overpotential(j, j_lim_anode) for j in j_values]
eta_conc_cathode = [concentration_overpotential(j, j_lim_cathode) for j in j_values]

# Ohmic losses
R_area = 0.1 * 1e-4  # Ω·m²
eta_ohm = [j * R_area for j in j_values]

# Compute operating voltages
V_operating = [E_Nernst - (eta_act_anode[i] + eta_act_cathode[i] + eta_conc_anode[i] + eta_conc_cathode[i] + eta_ohm[i]) for i in 1:length(j_values)]

# Plot polarization curve
plot(j_values ./ 1e4, V_operating, xlabel="Current Density (A/cm²)", ylabel="Cell Voltage (V)", title="Polarization Curve", legend=false, linewidth=2)

# Display results
println("Operating voltages for different current densities:")
for (i, j) in enumerate(j_values)
    println("j = ", j / 1e4, " A/cm² -> V = ", round(V_operating[i], digits=4), " V")
end
