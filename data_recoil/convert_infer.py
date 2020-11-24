#!/usr/bin/env python

import numpy as np
import root_numpy as rnp
import sys
from re import sub

if __name__ == '__main__':
    assert len(sys.argv) == 2

    froot = sys.argv[1]
    fnpy = sub('\.root$', '.npz', froot)
    print(froot)

    branches = ['pt', 'eta', 'phi', 'e', 'puppi', 'pdgid', 'hardfrac', 'cluster_idx', 'vtxid',
                'cluster_r', 'cluster_hardch_pt', 'cluster_puch_pt', 'npv','isolep']
    genz_branches = ['genZpt', 'genZphi']
    recz_branches = ['recZpt', 'recZphi']
    met_branches = ['genmet', 'genmetphi']
    recoil_branches = ['genUmag', 'genUphi']
    jet1_branches = ['genjet1%s'%b for b in ['pt', 'eta', 'phi', 'e']]
    jet2_branches = ['genjet2%s'%b for b in ['pt', 'eta', 'phi', 'e']]

    arr = rnp.root2array(filenames=[froot], treename='events', branches=branches+met_branches+recoil_branches+jet1_branches+jet2_branches+genz_branches+recz_branches)
    genz = np.stack([arr[k] for k in genz_branches], axis=-1)
    recz = np.stack([arr[k] for k in recz_branches], axis=-1)
    met = np.stack([arr[k] for k in met_branches], axis=-1)
    recoil = np.stack([arr[k] for k in recoil_branches], axis=-1)
    jet1 = np.stack([arr[k] for k in jet1_branches], axis=-1)
    jet2 = np.stack([arr[k] for k in jet2_branches], axis=-1)
    arrs = [np.stack(arr[k], axis=0) for k in branches]
    arr = np.stack(arrs, axis=-1)

    arr = arr.astype(np.float32)
    genz = genz.astype(np.float32)
    recz = recz.astype(np.float32)
    met = met.astype(np.float32)
    recoil = recoil.astype(np.float32)
    jet1 = jet1.astype(np.float32)
    jet2 = jet2.astype(np.float32)

    np.savez(fnpy, x=arr, met=met, recoil=recoil, jet1=jet1, jet2=jet2, genz=genz, recz=recz)

