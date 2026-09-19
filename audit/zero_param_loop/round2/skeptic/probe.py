import sys, json, numpy as np
sys.path.insert(0,"/mnt/disks/disk-socrateai-local-1/dualscale-wt-loop"); sys.path.insert(0,"/mnt/disks/disk-socrateai-local-1/dualscale-wt-loop/scripts")
sys.path.insert(0,"/tmp/claude-1501372770/-home-callensxavier-gmail-com-SocrateAI-Scientific-DualScaleSimulator/a29a4e9b-897f-4195-899e-9c4c27a67e9d/scratchpad/de")
import decisive_experiment as de, param_loop_sim as pls
np.random.seed(42)
def chi2_of_point(res):
    d=res["dark_energy"]
    zg=np.array(d["z_grid"]); DM=np.array(d["D_M_times_H0_grid"]); Hz=np.array(d["H_of_z_over_H0_grid"])
    xM={z:float(np.interp(z,zg,DM)) for z in de.DESI_Z_UNIQUE}; xH={z:float(1/np.interp(z,zg,Hz)) for z in de.DESI_Z_UNIQUE}
    base=[]
    for z,_,q in de.DESI_ROWS:
        base.append(xM[z] if q=="DM_over_rs" else xH[z] if q=="DH_over_rs" else (z*xM[z]**2*xH[z])**(1/3))
    bao=de.chi2_bao_given_base(np.array(base))
    DLq=(1+de.SN_Z)*np.interp(de.SN_Z,zg,DM); sn=de.chi2_sn_given_mu(5*np.log10(np.clip(DLq,1e-12,None)))
    return bao["chi2"], sn["chi2"], bao["chi2"]+sn["chi2"]
nom=dict(pls.DEFAULT_PARAMS)
probes={"nominal":{}}
for p in ["a_pot","b_pot","lambda_sym","mu_sym"]:
    probes[p+"_x10"]={p:nom[p]*10}; probes[p+"_x0.1"]={p:nom[p]/10}
out={}
for name,ov in probes.items():
    prm=dict(nom); prm.update(ov)
    r=pls.evaluate_point(prm)
    d=r["dark_energy"]; s=r["screening"]; pta=r["pta"]
    try: c=chi2_of_point(r)
    except Exception as e: c=("ERR",str(e))
    out[name]={"params":ov,"w0_cpl":d.get("w0_cpl_latetime"),"wa_cpl":d.get("wa_cpl_latetime"),"HzH0_z1":float(np.interp(1.0,d["z_grid"],d["H_of_z_over_H0_grid"])) if d.get("H_of_z_over_H0_grid") else None,
     "ssf":s["screening_suppression_factor"],"pcr":s["phi_center_ratio"],"stable":s["numerically_stable"],"pta_maxdev":pta["max_deviation_from_hd"],"chi2_bao_sn_total":c}
    print(name, json.dumps(out[name], default=str))
json.dump(out,open("probe_report.json","w"),indent=2,default=str)
