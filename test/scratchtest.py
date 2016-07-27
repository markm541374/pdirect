import pyximport
import numpy as np
pyximport.install(setup_args={ "include_dirs":[np.get_include()]})
import scratch
print scratch.vv()