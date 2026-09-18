"""Active-phase curvature and variable-phase, exactly conservative Gibbs splitting.
Diagnostic only: a finite local/multistart search is not proof of global stability.
"""
from __future__ import annotations
import argparse, ctypes, json, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, least_squares
from scipy.special import softmax

class CPA:
    def __init__(self, library: Path, temperature: float, pressure: float):
        self.lib=ctypes.CDLL(str(library.resolve())); self.t=temperature; self.p=pressure
        ptr=np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
        self.lib.split_mu.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.c_int,ptr,ptr,ptr]
        self.lib.split_original_flash.argtypes=[ctypes.c_double,ctypes.c_double,ptr,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int)]
        self.lib.split_last_error.restype=ctypes.c_char_p
        self.lib.split_feed.argtypes=[ptr]
        self.z=np.empty(5); self.lib.split_feed(self.z)
    def mu(self, x: np.ndarray, root: int):
        x=np.ascontiguousarray(x,dtype=float); m=np.empty(5); zz=np.empty(1)
        if self.lib.split_mu(self.t,self.p,int(root),x,m,zz):
            raise RuntimeError(self.lib.split_last_error().decode())
        return m,float(zz[0])
    def initial(self):
        n=np.zeros(15); r=(ctypes.c_int*3)(); s=ctypes.c_int()
        q=self.lib.split_original_flash(self.t,self.p,n,r,ctypes.byref(s))
        if q<1: raise RuntimeError(self.lib.split_last_error().decode())
        return n[:5*q].reshape(q,5),np.array(list(r)[:q]),bool(s.value)

def curvature(eos,x,root,h=1e-5):
    """Constrained Hessian in x[0:n-1], x[-1]=1-sum(x[:-1])."""
    d=len(x)-1; H=np.empty((d,d))
    for j in range(d):
        step=h*min(x[j],x[-1]); xp=x.copy(); xm=x.copy()
        xp[j]+=step; xp[-1]-=step; xm[j]-=step; xm[-1]+=step
        a,_=eos.mu(xp,root); b,_=eos.mu(xm,root)
        H[:,j]=((a[:-1]-a[-1])-(b[:-1]-b[-1]))/(2*step)
    asym=float(np.max(np.abs(H-H.T)))
    vals,vec=np.linalg.eigh((H+H.T)/2)
    v=np.r_[vec[:,0],-sum(vec[:,0])]; v/=np.linalg.norm(v)
    return float(vals[0]),v,asym

def state(eos,n,roots):
    if np.any(n<=0) or not np.all(np.isfinite(n)): raise ValueError('Positive finite inventory required')
    beta=n.sum(axis=1); x=n/beta[:,None]; mus=[]; zz=[]; curv=[]; errors=[]; directions=[]
    for xx,r in zip(x,roots):
        m,Z=eos.mu(xx,r); mus.append(m); zz.append(Z)
        a,v,_=curvature(eos,xx,r); b,_,_=curvature(eos,xx,r,5e-6)
        curv.append(b); errors.append(abs(a-b)); directions.append(v)
    mus=np.array(mus); curv=np.array(curv); errors=np.array(errors)
    uncertainty=np.maximum(1e-6,10*errors)
    # Negative curvature is a definite rejection. Near-zero curvature is UNKNOWN.
    local='UNSTABLE' if np.any(curv < -uncertainty) else ('LOCALLY_STABLE' if np.all(curv>uncertainty) else 'NEAR_ZERO_UNRESOLVED')
    return dict(n=n,x=x,beta=beta,roots=roots,mu=mus,Z=np.array(zz),g=float(np.sum(n*mus)),
                curvature=curv,curvature_error=errors,directions=np.array(directions),local_status=local,
                mass_error=float(np.max(np.abs(n.sum(axis=0)-eos.z))),
                fugacity_spread=float(np.max(np.ptp(mus,axis=0))))

def conservative_split(eos,n,roots,phase,v):
    """Symmetric daughter seed, exact component inventory; accept only G descent."""
    x=n[phase]/sum(n[phase]); amount=sum(n[phase]); base=state(eos,n,roots)['g']
    limit=min(x[abs(v)>1e-16]/abs(v[abs(v)>1e-16]))
    best=None
    for frac in (.05,.15,.3,.6):
        dx=frac*limit*v
        a=amount*.5*(x+dx); b=n[phase]-a
        trial=np.vstack([n[:phase],a,b,n[phase+1:]])
        rr=np.r_[roots[:phase],roots[phase],roots[phase],roots[phase+1:]]
        g=sum(float(nn @ eos.mu(nn/sum(nn),r)[0]) for nn,r in zip(trial,rr))
        if g < base-1e-13 and (best is None or g<best[2]):best=(trial,rr,g)
    if best is None: raise RuntimeError('Negative curvature did not produce a resolved conservative descent seed')
    return best

def minimize_allocations(eos,n,roots):
    """q unrestricted phase records; softmax of per-component allocations preserves z.
    Roots are attached to records, NOT unique O/G/W slots. q may exceed three.
    """
    q,nc=n.shape; a=n/eos.z; a=a/a.sum(axis=0)
    u0=np.log(np.maximum(a[:-1],1e-150)/np.maximum(a[-1],1e-150)).ravel()
    def eval_u(u):
        alpha=softmax(np.vstack([u.reshape(q-1,nc),np.zeros(nc)]),axis=0)
        nn=alpha*eos.z; xx=nn/nn.sum(axis=1)[:,None]
        mm=np.array([eos.mu(x,int(r))[0] for x,r in zip(xx,roots)])
        g=float(np.sum(nn*mm)); avg=np.sum(alpha*mm,axis=0)
        grad=(nn*(mm-avg))[:-1].ravel()
        return g,grad
    result=minimize(eval_u,u0,jac=True,method='BFGS',options={'gtol':2e-12,'maxiter':800})
    # Inventory-weighted gradients can hide trace-species chemical-potential
    # mismatch. Correct ALL phase/component mu equalities at exact fixed mass.
    def chemical_residual(u):
        aa=softmax(np.vstack([u.reshape(q-1,nc),np.zeros(nc)]),axis=0)
        mmass=aa*eos.z; xx=mmass/mmass.sum(axis=1)[:,None]
        mu=np.array([eos.mu(x,int(r))[0] for x,r in zip(xx,roots)])
        return (mu[:-1]-mu[-1]).ravel()
    corrected=least_squares(chemical_residual,result.x,jac='3-point',diff_step=1e-5,
        xtol=2e-13,ftol=2e-13,gtol=2e-13,max_nfev=180)
    u=corrected.x
    g_after,_=eval_u(u)
    if g_after > result.fun+1e-11: raise RuntimeError('Stationarity correction increased G')
    alpha=softmax(np.vstack([u.reshape(q-1,nc),np.zeros(nc)]),axis=0)
    nn=alpha*eos.z
    return nn,{'success':bool(result.success),'message':str(result.message),'nit':result.nit,
               'gradient_max':float(np.max(np.abs(result.jac))),'objective':float(result.fun),
               'chemical_corrector_success':bool(corrected.success),
               'chemical_residual_max':float(np.max(np.abs(corrected.fun))),
               'chemical_corrector_nfev':corrected.nfev}

def tpd_screen(eos,s):
    """Test BOTH candidate roots, including roots already active; no free-slot filter.
    Every active phase gets its own reference. Finite multistart is a screen only.
    """
    records=[]; ncomp=len(eos.z)
    for parent,(ref,xparent) in enumerate(zip(s['mu'],s['x'])):
        best=(0.,None,None); failures=0
        seeds=list(s['x'])+[eos.z,np.ones(ncomp)/ncomp]
        seeds += [np.eye(ncomp)[i]*.9+np.ones(ncomp)*.1/ncomp for i in range(ncomp)]
        for root in (0,1):
            def fun(u):
                x=softmax(np.r_[u,0.]); mu,_=eos.mu(x,root)
                dmu=mu-ref; f=float(x @ dmu)
                return f,(x*(dmu-f))[:-1]
            for x0 in seeds:
                u=np.log(np.maximum(x0[:-1],1e-100)/max(x0[-1],1e-100))
                try:
                    r=minimize(fun,u,jac=True,method='BFGS',options={'gtol':1e-9,'maxiter':100})
                    # A negative evaluated TPD is a witness even if optimizer stops early.
                    if not r.success: failures+=1
                    if np.isfinite(r.fun) and r.fun<best[0]: best=(float(r.fun),root,softmax(np.r_[r.x,0.]).tolist())
                except (RuntimeError,FloatingPointError): failures+=1
        records.append({'parent':parent,'minimum_tpd':best[0],'candidate_root':best[1],
                        'candidate_composition':best[2],'incomplete_searches':failures})
    return records

def serial(s):
    return {k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in s.items() if k!='directions'}

def run_case(lib,t,p,out,max_phases=6):
    started=time.monotonic(); eos=CPA(lib,t,p); n,roots,old_stable=eos.initial(); initial=state(eos,n,roots)
    report={'T_K':t,'P_MPa':p,'legacy_missing_phase_stable':old_stable,'initial':serial(initial),'splits':[]}
    s=initial
    for _ in range(max_phases-len(n)):
        if s['local_status']=='UNSTABLE':
            j=int(np.argmin(s['curvature'])); trigger='NEGATIVE_ACTIVE_CURVATURE'
            trial,rr,seed_g=conservative_split(eos,n,roots,j,s['directions'][j])
        else:
            witnesses=tpd_screen(eos,s)
            witness=min(witnesses,key=lambda r:r['minimum_tpd'])
            if witness['minimum_tpd']>=-1e-7: break
            j=witness['parent']; trigger='FINITE_DISTANCE_TPD'
            xc=np.array(witness['candidate_composition']); rc=witness['candidate_root']
            # Remove candidate inventory from its own parent, not from an empty slot.
            limit=min(n[j]/np.maximum(xc,1e-300)); best=None
            for fraction in (.001,.01,.1,.3):
                child=fraction*limit*xc; rest=n.copy();rest[j]-=child
                nn=np.vstack([rest,child]);rrr=np.r_[roots,rc]
                gg=sum(float(v @ eos.mu(v/sum(v),r)[0]) for v,r in zip(nn,rrr))
                if gg<s['g']-1e-13 and (best is None or gg<best[2]):best=(nn,rrr,gg)
            if best is None:
                report['split_stop']='NEGATIVE_TPD_NO_RESOLVED_DESCENT_SEED';break
            trial,rr,seed_g=best
        optimized,info=minimize_allocations(eos,trial,rr)
        nxt=state(eos,optimized,rr)
        if nxt['g'] > seed_g+1e-12 or nxt['mass_error']>1e-10:raise RuntimeError('Gibbs or mass invariant violated')
        report['splits'].append({'parent':j,'trigger':trigger,'seed_delta_G':seed_g-s['g'],
            'relaxed_delta_G':nxt['g']-s['g'],'optimizer':info,'state':serial(nxt)})
        n,roots,s=optimized,rr,nxt
    report['final']=serial(s)
    report['tpd']=tpd_screen(eos,s)
    has_negative=any(r['minimum_tpd'] < -1e-7 for r in report['tpd'])
    searches_incomplete=any(r['incomplete_searches'] for r in report['tpd'])
    report['equilibrium_status']=('REJECTED_UNSTABLE' if s['local_status']=='UNSTABLE' or has_negative else
        'UNRESOLVED_STATIONARITY' if s['fugacity_spread']>1e-6 else
        'UNRESOLVED_TPD_SEARCH' if searches_incomplete else 'SCREENED_NOT_GLOBALLY_PROVEN')
    report['FLOW_GATE']='BLOCKED';report['FIGURE7_REPRODUCTION']='BLOCKED'
    report['seconds']=time.monotonic()-started
    out.write_text(json.dumps(report,indent=2,allow_nan=False))
    print(t,p,'q=',len(n),'G decrease=',s['g']-initial['g'],'local=',s['local_status'],
          'mu=',s['fugacity_spread'],'tpd=',min([r['minimum_tpd'] for r in report['tpd']] or [0]),flush=True)
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('library',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--temperature',type=float,required=True);parser.add_argument('--pressure',type=float,required=True)
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    run_case(args.library,args.temperature,args.pressure,args.output)
