"""
Script to execute the entire pipeline of the dvm model
"""
import os
from pathlib import Path
import time
path_to_scripts = Path('~/model_dvm_carbon/Scripts/').expanduser()

start = time.time()

# Migration model with visual predation
VP_model = 0
if VP_model == 1:
    script = path_to_scripts / "dvm_model_vp.py"
    os.system("python " + str(script))
    script = path_to_scripts / "Resultss.py"
    os.system("python " + str(script))

# Migration model with annual variation of K
K_model = 0
if K_model == 1:
    script = path_to_scripts / "dvm_model_seasonal_v1.py"
    os.system("python " + str(script))
    script = path_to_scripts / "Results_env.py"
    os.system("python " + str(script))

# Migration model with annual variation of T and rho(T)
T_model = 1
if T_model == 1:
    script = path_to_scripts / "dvm_model_seasonal_v2.py"
    os.system("python " + str(script))
    script = path_to_scripts / "Results_envT.py"
    os.system("python " + str(script))

# Migration model with annual variation of irradiance
Ir_model = 0
if Ir_model == 1:
    script = path_to_scripts / "dvm_model_seasonal_v3.py"
    os.system("python " + str(script))
    script = path_to_scripts / "Results_envI.py"
    os.system("python " + str(script))

# Migration model with annual variation of K, T and irradiance
KTI_model = 0
if KTI_model == 1:
    script = path_to_scripts / "dvm_model_seasonal_v4.py"
    os.system("python " + str(script))
    script = path_to_scripts / "Results_envKTI.py"
    os.system("python " + str(script))

# Test the influence of size and taxonomy
size_descriptors = 0
if size_descriptors == 1:
    script = path_to_scripts / "size_descriptors.py"
    os.system("python " + str(script))

end = time.time()
print("Execution took", (end-start)/60, "min")