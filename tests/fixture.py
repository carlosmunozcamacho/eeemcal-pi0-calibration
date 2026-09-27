# SPDX-License-Identifier: LGPL-3.0-or-later
"""Small deterministic truth-free reconstructed fixture, generated on demand."""
from pathlib import Path
import numpy as np
import awkward as ak
import uproot


def reco(path,n=2200,missing=None):
    rng=np.random.default_rng(741);cells=np.arange(8,dtype=np.uint64)+2**54+1
    m=rng.normal(.1285,.0055,n);bg=rng.random(n)<.22;m[bg]=rng.uniform(.04,.24,bg.sum())
    e=rng.uniform(.7,3.,(n,2));theta=np.arccos(1-m*m/(2*np.prod(e,axis=1)))
    pos=np.zeros((n,2,3));pos[:,:,2]=-1850;pos[:,1,0]=1850*np.tan(theta)
    hitix=[rng.choice(8,4,replace=False).tolist() for _ in range(n)]
    data={}
    def put(c,k,v,dtype=None): data[c+'/'+c+'.'+k]=ak.Array(np.asarray(v,dtype=dtype)) if dtype else ak.Array(v)
    cl='EcalEndcapNClusters';rh='EcalEndcapNRecHits'
    put(cl,'energy',e)
    for q,d in enumerate('xyz'): put(cl,'position.'+d,pos[:,:,q])
    for k,v in [('nhits',[2,2]),('hits_begin',[0,2]),('hits_end',[2,4]),('hitContributions_begin',[0,2]),('hitContributions_end',[2,4])]:put(cl,k,np.tile(v,(n,1)))
    contrib=np.column_stack([.7*e[:,0],.3*e[:,0],.7*e[:,1],.3*e[:,1]])
    data['_'+cl+'_hitContributions']=ak.Array(contrib)
    put('_'+cl+'_hits','index',np.tile([0,1,2,3],(n,1)));put('_'+cl+'_hits','collectionID',np.full((n,4),2))
    put(rh,'cellID',cells[np.array(hitix)],np.uint64);put(rh,'energy',contrib)
    for q,d in enumerate('xyz'):
        values=np.array(hitix)*20 if q==0 else np.zeros((n,4)) if q==1 else np.full((n,4),-1850)
        put(rh,'position.'+d,values)
    put('EventHeader','eventNumber',np.arange(n).reshape(-1,1));put('EventHeader','runNumber',np.ones((n,1),dtype=int))
    empty=ak.unflatten(np.array([],dtype=np.int32),np.zeros(n,dtype=int))
    for rel in ['cluster','track']:
        for member in ['index','collectionID']: put('_EcalEndcapNTrackClusterMatches_'+rel,member,empty)
    if missing: data.pop(missing)
    with uproot.recreate(path) as f:
        f.mktree('events',data)
        f['podio_metadata']=ak.Array([{'events___CollectionTypeInfo':[{'name':name,'collectionID':i+1} for i,name in enumerate(['EcalEndcapNClusters','EcalEndcapNRecHits','CentralCKFTracks'])]}])
    return cells

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--events',type=int,default=2200)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True);reco(a.output/'toy.root',a.events)
    (a.output/'manifest.txt').write_text(str((a.output/'toy.root').resolve())+'\n')
