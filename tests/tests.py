#!/usr/bin/env python
# Path hack.
import sys, os
sys.path.insert(0, os.path.abspath('..'))

print(sys.version_info)

if sys.version_info < (3, 0):
    print("Using Python 2.7")
    # Just test whether the import works
    import kimg_pipe
else:
    print("Using Python 3.5")
    from kimg_pipe import kimg_pipe

patient = kimg_pipe.freeCoG(subj = 'S1', hem = 'lh', fs_dir='', subj_dir='.')

# Test whether mayavi works
import mayavi
from mayavi import mlab

mlab.test_mesh()

# Test import of plotting
import kimg_pipe.plotting as plotting
import kimg_pipe.SupplementalScripts
