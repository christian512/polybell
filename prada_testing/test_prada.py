from linearbell.utils import get_deterministic_behaviors
from linearbell.panda_helper import write_known_vertices, read_inequalities
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import time

input_output_configs = [(4,3,2,2)]
number_of_inequalities_to_find = [12480]
recursion_levels_to_test = [0]
for config, num_ineqs in zip(input_output_configs, number_of_inequalities_to_find):
    times = []
    for recursion_level in recursion_levels_to_test:
        start_time = time.time()
        # same number of outputs
        assert config[2] == config[3]
        vertices = get_deterministic_behaviors(range(config[0]), range(config[1]), range(config[2]))
        relabels_dets = np.loadtxt('data/relabels_dets/{}{}{}{}.gz'.format(config[0], config[1], config[2], config[3])).astype(int)
        write_known_vertices(vertices, 'randa_test_file.ext', relabellings_vertices=relabels_dets)
        # run randa
        cmd = 'randa randa_test_file.ext -t 1 -r ' + str(recursion_level) + ' > out.ine'
        out = subprocess.run(cmd, shell=True)
        inequalities = read_inequalities('out.ine').astype(float)
        times.append(time.time() - start_time)
        assert inequalities.shape[0] == num_ineqs, inequalities.shape[0]
    plt.scatter(recursion_levels_to_test, times, label=str(config))
    
plt.title('Speed comparison Recursive Panda w/o recalculation checks ')
plt.xlabel('Recursion depth')
plt.ylabel('Time (s)')
plt.legend()
plt.show()
