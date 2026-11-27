import math, numpy as np
from tss.l3_signal.novel import drift_crb as D
from tss.l3_signal.novel.route_one_estimator import estimate_route_one
from tss.l3_signal.novel.route_one_joint import estimate_route_one_joint, joint_model
from tss.constants import C_LIGHT
F,P,MU,W,VS=4.988e6,529.8e3,0.9537,0.926,7600.0
T=(np.arange(209)-104)/45.0
S=2*F*VS**2/(C_LIGHT*P); K=MU*P/VS**2; v=25.0
for tau in (0.5,2.0):
  R=np.exp(-np.abs(T[:,None]-T[None,:])/tau); L=np.linalg.cholesky(R+1e-12*np.eye(T.size))
  for snr in (20.0,10.0):
    crb=D.route_one_crb(D.fisher_route_one(T,snr_db=snr,width_s=W,f_hz=F,mu_sat=MU,p_prime_m=P,v_s_mps=VS,drift_mps=v,tau_c_s=tau)).sigma["v"]
    rng=np.random.default_rng(2); A=10**(snr/20); ep,ej,f=[],[],0
    for i in range(200):
      g=L@((rng.standard_normal(T.size)+1j*rng.standard_normal(T.size))/math.sqrt(2))
      y=A*joint_model(T,t0_s=0,drift_mps=v,width_s=W,slope_hz_per_s=S,kappa_s2_per_m=K)*g+(rng.standard_normal(T.size)+1j*rng.standard_normal(T.size))/math.sqrt(2)
      e=estimate_route_one(T,y,f_hz=F,group_path_m=P,mu_sat=MU,v_s_mps=VS)
      if not e.separable: continue
      f+=1; j=estimate_route_one_joint(T,y,f_hz=F,group_path_m=P,mu_sat=MU,v_s_mps=VS,t0_start_s=e.t_power_peak_s,drift_start_mps=e.drift_mps)
      ep.append(e.drift_mps-v); ej.append(j.drift_mps-v if j.ok else np.nan)
    ep,ej=np.array(ep),np.array(ej); r=lambda x: math.sqrt(np.nanmean(x**2))
    print(f"tau {tau} {snr:.0f} dB formed {f}/200 fading CRB {crb:.2f} pinned {r(ep):.2f} ({r(ep)/crb:.2f}x) joint {r(ej):.2f} ({r(ej)/crb:.2f}x) robust-sd joint {1.4826*np.nanmedian(np.abs(ej-np.nanmedian(ej))):.2f}")
