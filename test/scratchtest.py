import pyximport
import numpy as np
pyximport.install(setup_args={ "include_dirs":[np.get_include()]})
import scratch
print scratch.foo(4,3)
scratch.bar([1,23])