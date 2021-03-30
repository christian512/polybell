from linearbell.utils import get_deterministic_behaviors
from linearbell.panda_helper import write_known_vertices, read_inequalities
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import time

input_output_configs = [(2, 2, 2, 2), (3,3,2,2), (4,3,2,2)]
times_no_symm = []
times_symm = []

for config in input_output_configs:
    start_time = time.time()
    # same number of outputs
    assert config[2] == config[3]
    vertices = get_deterministic_behaviors(range(config[0]), range(config[1]), range(config[2]))
    write_known_vertices(vertices, 'randa_test_file.ext')
    # run randa
    cmd = 'randa randa_test_file.ext -t 1 -r ' + str(0) + ' > out.ine'
    out = subprocess.run(cmd, shell=True)
    times_no_symm.append(time.time() - start_time)
plt.scatter(list(range(len(times_no_symm))), times_no_symm, label="No symmetry check on top-level")

for config in input_output_configs:
    start_time = time.time()
    # same number of outputs
    assert config[2] == config[3]
    vertices = get_deterministic_behaviors(range(config[0]), range(config[1]), range(config[2]))
    relabels_dets = np.loadtxt('../data/relabels_dets/{}{}{}{}.gz'.format(config[0], config[1], config[2], config[3])).astype(int)
    write_known_vertices(vertices, 'randa_test_file.ext', relabellings_vertices=relabels_dets)
    # run randa
    cmd = 'randa randa_test_file.ext -t 1 -r ' + str(0) + ' > out.ine'
    out = subprocess.run(cmd, shell=True)
    times_symm.append(time.time() - start_time)
plt.scatter(list(range(len(times_symm))), times_symm, label="Using Symmetry check on top-level")
    
plt.title('Speed comparison Adjacency Decomposition with and w/o equivalence checks')
plt.xlabel('Different configurations (increaing dimensions)')
plt.ylabel('Time (s)')
plt.legend()
plt.show()
