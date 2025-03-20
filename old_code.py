import cantera as ct
from scipy.optimize import differential_evolution, minimize
import os

def chemical_heat_release_H2(gas, reactant_massfracs, temperature, pressure):
    """
    Compute the heat of combustion for H2

            H2 + 0.5*O2 --> H2O

    Return in J/kg
    """
    # Store here all the chemical info needed
    hf_N2 = 0.0                     # J/mol
    hf_H = 218.00  * (10**3)        # J/mol
    hf_O2 = 0.0                     # J/mol
    hf_O = 249.18 * (10**3)         # J/mol
    hf_OH = 38.99 * (10**3)         # J/mol
    hf_H2 = 0.0                     # J/mol
    hf_H2O = -241.83 * (10**3)      # J/mol
    hf_HO2 = 2.09 * (10**3)         # J/mol
    hf_H2O2 = -136.11 * (10**3)     # J/mol

    MW_N2 = gas.molecular_weights[0] * (10**-3)      # kg/mol
    MW_H = gas.molecular_weights[1] * (10**-3)       # kg/mol
    MW_O2 = gas.molecular_weights[2] * (10**-3)      # kg/mol
    MW_O = gas.molecular_weights[3] * (10**-3)       # kg/mol
    MW_OH = gas.molecular_weights[4] * (10**-3)      # kg/mol
    MW_H2 = gas.molecular_weights[5] * (10**-3)      # kg/mol
    MW_H2O = gas.molecular_weights[6] * (10**-3)     # kg/mol
    MW_HO2 = gas.molecular_weights[7] * (10**-3)     # kg/mol
    MW_H2O2 = gas.molecular_weights[8] * (10**-3)    # kg/mol

    hf_N2_mass = hf_N2 / MW_N2         # J/kg
    hf_H_mass = hf_H / MW_H            # J/kg
    hf_O2_mass = hf_O2 / MW_O2         # J/kg
    hf_O_mass = hf_O / MW_O            # J/kg
    hf_OH_mass = hf_OH / MW_OH         # J/kg
    hf_H2_mass = hf_H2 / MW_H2         # J/kg
    hf_H2O_mass = hf_H2O / MW_H2O      # J/kg
    hf_HO2_mass = hf_HO2 / MW_HO2      # J/kg
    hf_H2O2_mass = hf_H2O2 / MW_H2O2   # J/kg

    # Get reactant mass fracs
    yN2_unburned = reactant_massfracs[0]
    yH_unburned = reactant_massfracs[1]
    yO2_unburned = reactant_massfracs[2]
    yO_unburned = reactant_massfracs[3]
    yOH_unburned = reactant_massfracs[4]
    yH2_unburned = reactant_massfracs[5]
    yH2O_unburned = reactant_massfracs[6]
    yHO2_unburned = reactant_massfracs[7]
    yH2O2_unburned = reactant_massfracs[8]

    # Get equilibrium state here
    gas.TP = temperature, pressure
    gas.equilibrate('TV')
    equil_massfracs = gas.Y
    yN2_burned = equil_massfracs[0]
    yH_burned = equil_massfracs[1]
    yO2_burned = equil_massfracs[2]
    yO_burned = equil_massfracs[3]
    yOH_burned = equil_massfracs[4]
    yH2_burned = equil_massfracs[5]
    yH2O_burned = equil_massfracs[6]
    yHO2_burned = equil_massfracs[7]
    yH2O2_burned = equil_massfracs[8]

    # Compute the heats of formation of reactant and product states (J/kg)
    hf1 = yN2_unburned*hf_N2_mass + yH_unburned*hf_H_mass + yO2_unburned*hf_O2_mass + yO_unburned*hf_O_mass \
        + yOH_unburned*hf_OH_mass + yH2_unburned*hf_H2_mass + yH2O_unburned*hf_H2O_mass + yHO2_unburned*hf_HO2_mass \
        + yH2O2_unburned*hf_H2O2_mass

    hf2 = yN2_burned*hf_N2_mass + yH_burned*hf_H_mass + yO2_burned*hf_O2_mass + yO_burned*hf_O_mass \
        + yOH_burned*hf_OH_mass + yH2_burned*hf_H2_mass + yH2O_burned*hf_H2O_mass + yHO2_burned*hf_HO2_mass \
        + yH2O2_burned*hf_H2O2_mass
    qc = hf2 - hf1
    return qc


def compute_CJ_solution(reactant_density, reactant_pressure, reactant_gamma, qc):
    """
    Compute the unburned Mach number required
    for the upper CJ point
    """
    rho = reactant_density
    pres = reactant_pressure
    gamma = reactant_gamma
    q_hat_c = -1*qc * (rho/pres)

    def Mu_CJ(gamma, qhatc):
        return (1 + (((gamma**2 - 1)*qhatc)/(gamma))*(1 + (1 + (2*gamma)/((gamma**2 - 1)*qhatc))**(0.5)))**(0.5)

    def v_CJ(gamma, qhatc):
        return 1 + (qhatc*((gamma - 1)/(gamma)))*(1 - (1 + (2*gamma)/((gamma**2 - 1)*qhatc))**(0.5))

    def p_CJ(gamma, qhatc):
        return 1 + (qhatc*(gamma - 1))*(1 + (1 + (2*gamma)/((gamma**2 - 1)*qhatc))**(0.5))

    mach_plus_CJ = Mu_CJ(gamma, q_hat_c)
    v_hat_CJ = v_CJ(gamma, q_hat_c)
    p_hat_CJ = p_CJ(gamma, q_hat_c)

    return mach_plus_CJ, v_hat_CJ, p_hat_CJ


def main():
    # Define the mixture
    cwd = os.getcwd()
    if cwd == '/home/ltt/PeleC/Exec':
        gas = ct.Solution('/home/ltt/PeleC/Exec/Production/DetonationTest/HydrogenCK/hydrogen.yaml')
    elif cwd == '/home/ltt/PeleC/Exec/Production/DetonationTest':
        gas = ct.Solution('HydrogenCK/hydrogen.yaml')
    else:
        print("Need to be in correct starting directory")
        exit()

    # ---------------------------------
    # Set the pressure and temperature
    reac_pressure = 1.0 * ct.one_atm
    reac_temperature = 298.15  # Kelvin
    reac_massfracs = 'H2:0.1119, O2:0.88809'
    # ---------------------------------

    # Now set initial state
    gas.TPY = reac_temperature, reac_pressure, reac_massfracs
    reac_massfracs = gas.Y
    reac_cp = gas.cp_mass
    reac_cv = gas.cv_mass
    reac_gamma = reac_cp / reac_cv
    reac_rho = gas.density_mass
    reac_a = gas.sound_speed
    print("Specific heat ratio (gamma) for unburned gas: {:.2f}".format(reac_gamma))

    # Here goes ye olde optimization algorithm;
    # put the objective function here
    def full_CJ_algorithm(x):
        """
        The optimization algorithm will tweak the thermodynamic
        state of the burned state until convergence
        """
        bTemp, bPres = x
        # reac_temperature, reac_pressure, reac_massfracs, reac_gamma = args
        qc = chemical_heat_release_H2(gas, reac_massfracs, bTemp, bPres)
        temp_gamma = gas.cp_mass / gas.cv_mass
        mach_plus_CJ, v_hat_CJ, p_hat_CJ = compute_CJ_solution(reac_rho, reac_pressure, temp_gamma, qc)
        rho_burned = reac_rho / v_hat_CJ
        p_burned = p_hat_CJ * reac_pressure
        R = (ct.gas_constant) / (gas.mean_molecular_weight)
        temp_burned = p_burned / (rho_burned * R)
        return ((p_burned - bPres)**2)**(0.5) + ((temp_burned - bTemp)**2)**(0.5)

    # Run the algorithm
    result = minimize(full_CJ_algorithm, x0=[3600.0, 1_850_000.0], method='Powell', tol=1e-6, options= {'xtol':1e-3, 'maxfev':50000})
    if result.success == False:
        print("Optimization algorithm was not successful - exiting...")
        exit()
    final_temp, final_pres = result.x

    # Re-compute with proper thermodynamic state for the products
    qc_final = chemical_heat_release_H2(gas, reac_massfracs, final_temp, final_pres)
    print("--------------------------------")
    print("Heat release of hydrogen (H2) - oxygen (O2) combustion: {:.2f} MJ/kg".format(qc_final/(10**6)))
    print("--------------------------------")
    mach_plus_CJ, v_hat_CJ, p_hat_CJ = compute_CJ_solution(reac_rho, reac_pressure, reac_gamma, qc_final)
    print(f"Unburned Mach for CJ state is {mach_plus_CJ:.2f}")
    print(f"rho hat for CJ state is {1.0/v_hat_CJ:.2f}")
    print(f"temp hat for CJ state is {final_temp/reac_temperature:.2f}")
    print(f"p hat for CJ state is {p_hat_CJ:.2f}")
    final_a = gas.sound_speed

    # Compute extra info after solution
    print("--------------------------------")
    final_rho = reac_rho / v_hat_CJ
    print(f"rho burned for CJ state is {final_rho:.2f} kg/m^3")
    print(f"P burned for CJ state is {final_pres/reac_pressure:.2f} atm ({final_pres:.2f} pa)")
    print(f"Temp burned for CJ state is {final_temp:.1f} K")

    print("--------------------------------")
    final_massfracs = gas.Y
    print(f"yN2 burned for CJ state is {final_massfracs[0]:.5e}")
    print(f"yH burned for CJ state is {final_massfracs[1]:.5e}")
    print(f"yO2 burned for CJ state is {final_massfracs[2]:.5f}")
    print(f"yO burned for CJ state is {final_massfracs[3]:.5f}")
    print(f"yOH burned for CJ state is {final_massfracs[4]:.5f}")
    print(f"yH2 burned for CJ state is {final_massfracs[5]:.5f}")
    print(f"yH2O burned for CJ state is {final_massfracs[6]:.5f}")
    print(f"yHO2 burned for CJ state is {final_massfracs[7]:.5e}")
    print(f"yH2O2 burned for CJ state is {final_massfracs[8]:.5e}")

    print("--------------------------------")
    lab_unburned_vel = 0.0  # By definition here
    wave_unburned_vel = mach_plus_CJ * reac_a
    lab_det_speed = wave_unburned_vel
    wave_det_speed = 0.0    # Again, by definition of the frame
    wave_burned_vel = (1.0) * final_a
    lab_burned_vel = -1.0*wave_burned_vel + wave_unburned_vel
    lab_burned_mach_number = lab_burned_vel / final_a
    print("Wave Frame Velocities")
    print(f"Wave frame reactant velocity: {wave_unburned_vel:.1f} m/s")
    print(f"Wave frame detonation velocity: {wave_det_speed:.1f} m/s")
    print(f"Wave frame product velocity: {wave_burned_vel:.1f} m/s")
    print(f"Wave frame reactant Mach number: {mach_plus_CJ:.2f}")
    print(f"Wave frame detonation Mach number: {0.0:.2f}")
    print(f"Wave frame product Mach number: {1.0:.2f}")

    print("--------------------------------")
    print("Lab Frame Velocities")
    print(f"Lab frame reactant velocity: {lab_unburned_vel:.1f} m/s")
    print(f"Lab frame detonation velocity: {lab_det_speed:.1f} m/s")
    print(f"Lab frame product velocity: {lab_burned_vel:.1f} m/s")
    print(f"Lab frame reactant Mach number: {0.0:.2f}")
    print(f"Lab frame detonation Mach number: {mach_plus_CJ:.2f}")
    print(f"Lab frame product Mach number: {lab_burned_mach_number:.2f}\n")


if __name__ == "__main__":
    main()
