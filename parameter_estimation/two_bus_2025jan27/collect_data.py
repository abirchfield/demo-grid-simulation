import json
from esa import SAW
import numpy as np
import os
path = os.path.abspath(r"TwoBusJan27.PWB")
saw = SAW(path)
cases = []
for i in range(50):
    case = {}
    H = 2 + i*.1
    print(i,H)
    case["H"] = H
    genrou_frame = saw.GetParametersMultipleElement('MachineModel_GENROU', ["BusNum", "GenID", "TSH"])
    genrou_frame.at[1,"TSH"] = H
    saw.change_parameters_multiple_element_df("MachineModel_GENROU", genrou_frame)
    result = saw.RunScriptCommand(f'TSSolve("My Transient Contingency",[0,10,0.0041666666,NO]);')

    result = saw.TSGetContingencyResults("My Transient Contingency", ['Plot ''Main'''], 0, 10)
    rtab = result[1].to_numpy()
    case["time"] = list(rtab[:,0])
    case["headers"] = ["Vangle1", "Vangle2", "Vmag1", "Vmag2"]
    d = case["data"] = []
    for j in range(rtab.shape[1]):
        if j == 0: continue
        d.append(list(rtab[:,j]))
    cases.append(case)
with open("results.json", "w") as f:
    f.write(json.dumps(cases, indent=0))