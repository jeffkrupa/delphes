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
                'cluster_r', 'cluster_hardch_pt', 'cluster_puch_pt', 'npv']
    met_branches = ['genmet', 'genmetphi']
    jet1_branches = ['genjet1%s'%b for b in ['pt', 'eta', 'phi', 'e']]
    jet2_branches = ['genjet2%s'%b for b in ['pt', 'eta', 'phi', 'e']]

    arr = rnp.root2array(filenames=[froot], treename='events', branches=branches+met_branches+jet1_branches+jet2_branches)
    met = np.stack([arr[k] for k in met_branches], axis=-1)
    jet1 = np.stack([arr[k] for k in jet1_branches], axis=-1)
    jet2 = np.stack([arr[k] for k in jet2_branches], axis=-1)
    arrs = [np.stack(arr[k], axis=0) for k in branches]
    arr = np.stack(arrs, axis=-1)

    arr = arr.astype(np.float16)
    met = met.astype(np.float16)
    jet1 = jet1.astype(np.float16)
    jet2 = jet2.astype(np.float16)

    np.savez(fnpy, x=arr, met=met, jet1=jet1, jet2=jet2)
