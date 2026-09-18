"""Regression gates for the algorithm, not for physical phase-diagram validation."""
import json, sys
from pathlib import Path
import numpy as np
from audit import CPA, curvature, conservative_split, minimize_allocations, state, run_case, tpd_screen

class RegularSolution:
    z=np.array([.5,.5])
    def __init__(self,chi):self.chi=chi
    def mu(self,x,root):
        return np.log(x)+self.chi*np.array([x[1]**2,x[0]**2]),1.

lib=Path(sys.argv[1]); out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
for chi,expected in ((0,4.),(3,-2.)):
    value,_,_=curvature(RegularSolution(chi),np.array([.5,.5]),0)
    assert abs(value-expected)<1e-6,(chi,value)
# Binary regular-solution demixing: conservative seed and general phase allocation.
e=RegularSolution(3);n=np.array([[.5,.5]]);r=np.array([0]);s=state(e,n,r)
v=s['directions'][0]; seed,rr,g=conservative_split(e,n,r,0,v)
nn,info=minimize_allocations(e,seed,rr);q=state(e,nn,rr)
assert q['mass_error']<1e-14 and q['g']<s['g'] and q['fugacity_spread']<1e-8
assert q['local_status']=='LOCALLY_STABLE'
# A locally convex occupied liquid can still have a finite-distance TPD instability.
e=RegularSolution(3);e.z=np.array([.2,.8]);s=state(e,e.z[None,:],np.array([0]))
assert s['local_status']=='LOCALLY_STABLE'
assert tpd_screen(e,s)[0]['minimum_tpd'] < -1e-3
summary=[]
for t,p in ((588.,11.423873867918125),(590.,11.753904875999448),(593.,12.266165902108284)):
    report=run_case(lib,t,p,out/f'{t:g}.json')
    a,b=report['initial'],report['final'];count=len(b['beta'])
    assert report['legacy_missing_phase_stable']
    if t==588:
        assert a['local_status']=='LOCALLY_STABLE' and count==3
    else:
        assert a['local_status']=='UNSTABLE' and count>=4
        assert b['g']<a['g']-1e-10
        assert b['local_status']=='LOCALLY_STABLE'
        assert b['fugacity_spread']<1e-6 and b['mass_error']<1e-10
        # Reverse the unstable eigendirection: same two daughters in opposite order.
        eos=CPA(lib,t,p);n,roots,_=eos.initial();st=state(eos,n,roots)
        j=int(np.argmin(st['curvature']));nn,rr,_=conservative_split(eos,n,roots,j,-st['directions'][j])
        nn,_=minimize_allocations(eos,nn,rr);opposite=state(eos,nn,rr)
        assert abs(opposite['g']-b['g'])<1e-10
        xx=np.array(b['x']);yy=opposite['x']
        assert max(min(np.max(np.abs(v-w)) for w in yy) for v in xx)<1e-6
        report['reversed_direction_G_gap']=abs(opposite['g']-b['g'])
        (out/f'{t:g}.json').write_text(json.dumps(report,indent=2))
    summary.append({'T_K':t,'initial_phases':len(a['beta']),'final_phases':count,
        'initial_min_curvature':min(a['curvature']),'final_min_curvature':min(b['curvature']),
        'delta_G_over_nRT':b['g']-a['g'],'mass_error':b['mass_error'],
        'max_log_fugacity_spread':b['fugacity_spread'],'algorithm_regression':'PASS',
        'equilibrium_status':report['equilibrium_status'],'FORMAL_FLOW':'BLOCKED'})
(out/'summary.json').write_text(json.dumps(summary,indent=2))
print('ALGORITHM_REGRESSION_PASS; FIGURE7_REPRODUCTION=BLOCKED; FORMAL_FLOW=BLOCKED')
