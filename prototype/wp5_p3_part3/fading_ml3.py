import math, sys, time, numpy as np
from scipy.optimize import minimize
from scipy.linalg import cho_factor, cho_solve
from tss.l3_signal.novel import drift_crb as D
from tss.l3_signal.novel.route_one_estimator import estimate_route_one
from tss.l3_signal.novel.route_one_joint import joint_model
from tss.constants import C_LIGHT
F,P,MU,W,VS=4.988e6,529.8e3,0.9537,0.926,7600.0
T=(np.arange(209)-104)/45.0
S=2*F*VS**2/(C_LIGHT*P); K=MU*P/VS**2; v=25.0
LAG=np.abs(T[:,None]-T[None,:]); NOISE=1.0  # the tier measures it from the negative-beat half
def nll(q,y,tau_known=None):
    lp,t0,vv,lw=q[:4]; tau=tau_known if tau_known else math.exp(q[4])
    s=joint_model(T,t0_s=t0,drift_mps=vv,width_s=math.exp(lw),slope_hz_per_s=S,kappa_s2_per_m=K)
    C=math.exp(lp)*(s[:,None]*np.exp(-LAG/tau)*np.conj(s)[None,:])+NOISE*np.eye(T.size)
    try: cf=cho_factor(C,lower=True)
    except np.linalg.LinAlgError: return 1e30
    return float(2*np.sum(np.log(np.abs(np.diag(cf[0]))))+np.real(np.vdot(y,cho_solve(cf,y))))
def fit(y,e,tau_known=None):
    p=np.abs(y)**2
    tail=[] if tau_known else [math.log(0.5)]
    best=(math.inf,None)
    for lw in np.log([0.4,0.6,0.9,1.3,1.9]):
      s2=np.abs(joint_model(T,t0_s=e.t_power_peak_s,drift_mps=0,width_s=math.exp(lw),slope_hz_per_s=S,kappa_s2_per_m=K))**2
      lp=math.log(max(np.sum(p-NOISE)/np.sum(s2),1e-3))
      for t0 in np.linspace(-0.6,0.6,13):
        for vv in np.linspace(-100,150,26):
          q=[lp,t0,vv,lw]+tail; f=nll(q,y,tau_known)
          if f<best[0]: best=(f,q)
    r=minimize(nll,best[1],args=(y,tau_known),method="Nelder-Mead",options={"maxiter":4000,"xatol":1e-6,"fatol":1e-6})
    r=minimize(nll,r.x,args=(y,tau_known),method="Nelder-Mead",options={"maxiter":4000,"xatol":1e-7,"fatol":1e-7})
    return r.x[2], r
if __name__=="__main__":
  n=int(sys.argv[1]) if len(sys.argv)>1 else 30
  for tau in (0.5,2.0):
    Lc=np.linalg.cholesky(np.exp(-LAG/tau)+1e-12*np.eye(T.size))
    for snr in (20.0,10.0):
      crb=D.route_one_crb(D.fisher_route_one(T,snr_db=snr,width_s=W,f_hz=F,mu_sat=MU,p_prime_m=P,v_s_mps=VS,drift_mps=v,tau_c_s=tau)).sigma["v"]
      rng=np.random.default_rng(2); A=10**(snr/20); ep,ek,ef=[],[],[]; t1=time.time()
      for i in range(n):
        g=Lc@((rng.standard_normal(T.size)+1j*rng.standard_normal(T.size))/math.sqrt(2))
        y=A*joint_model(T,t0_s=0,drift_mps=v,width_s=W,slope_hz_per_s=S,kappa_s2_per_m=K)*g+(rng.standard_normal(T.size)+1j*rng.standard_normal(T.size))/math.sqrt(2)
        e=estimate_route_one(T,y,f_hz=F,group_path_m=P,mu_sat=MU,v_s_mps=VS)
        if not e.separable: continue
        ep.append(e.drift_mps-v); ek.append(fit(y,e,tau)[0]-v); ef.append(fit(y,e)[0]-v)
      r=lambda x: math.sqrt(np.mean(np.array(x)**2)); rs=lambda x: 1.4826*np.median(np.abs(np.array(x)-np.median(x)))
      print(f"tau {tau} {snr:.0f} dB formed {len(ep)}/{n} CRB {crb:.2f} | pinned {r(ep):.2f} | tau known {r(ek):.2f} ({r(ek)/crb:.2f}x, robust {rs(ek):.2f}) | tau free {r(ef):.2f} ({r(ef)/crb:.2f}x, robust {rs(ef):.2f}) [{time.time()-t1:.0f}s]",flush=True)
