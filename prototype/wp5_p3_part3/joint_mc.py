import sys, math, numpy as np
from tss.l3_signal.novel import drift_crb as D
from tss.l3_signal.novel.route_one_estimator import estimate_route_one
from tss.l3_signal.novel.route_one_joint import estimate_route_one_joint, joint_model
from tss.constants import C_LIGHT
F,P,MU,W,VS=4.988e6,529.8e3,0.9537,0.926,7600.0
T=(np.arange(209)-104)/45.0
S=2*F*VS**2/(C_LIGHT*P); K=MU*P/VS**2
v_true=25.0
n=int(sys.argv[1]) if len(sys.argv)>1 else 200
for snr in (20.0,15.0,10.0):
    crb=D.route_one_crb(D.fisher_route_one(T,snr_db=snr,width_s=W,f_hz=F,mu_sat=MU,p_prime_m=P,v_s_mps=VS,drift_mps=v_true)).sigma["v"]
    rng=np.random.default_rng(1)
    A=10**(snr/20); ep,ej,sj,formed,fail=[],[],[],0,0
    for i in range(n):
        y=A*joint_model(T,t0_s=0.0,drift_mps=v_true,width_s=W,slope_hz_per_s=S,kappa_s2_per_m=K)*np.exp(1j*rng.uniform(0,2*np.pi))
        y=y+(rng.standard_normal(T.size)+1j*rng.standard_normal(T.size))/math.sqrt(2)
        est=estimate_route_one(T,y,f_hz=F,group_path_m=P,mu_sat=MU,v_s_mps=VS)
        if not est.separable: continue
        formed+=1
        jt=estimate_route_one_joint(T,y,f_hz=F,group_path_m=P,mu_sat=MU,v_s_mps=VS,t0_start_s=est.t_power_peak_s,drift_start_mps=est.drift_mps)
        ep.append(est.drift_mps-v_true)
        if jt.ok: ej.append(jt.drift_mps-v_true); sj.append(jt.drift_sigma_mps)
        else: fail+=1
    ep,ej=np.array(ep),np.array(ej)
    rms=lambda e: math.sqrt(np.mean(e**2))
    print(f"{snr:4.0f} dB formed {formed}/{n} joint_fail {fail} CRB {crb:.3f} | pinned rms {rms(ep):.3f} ({rms(ep)/crb:.2f}x) bias {ep.mean():+.3f} | joint rms {rms(ej):.3f} ({rms(ej)/crb:.2f}x) bias {ej.mean():+.3f} mean sigma {np.mean(sj):.3f}")
