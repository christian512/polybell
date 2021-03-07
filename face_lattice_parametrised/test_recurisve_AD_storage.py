from linearbell.utils import get_deterministic_behaviors, get_configs, get_parametrisation_configs, parametrise_behavior
import numpy as np
from polytope import ParamPolytope

inputs_a = range(4)
inputs_b = range(3)
outputs = range(2)
configs = get_configs(inputs_a, inputs_b, outputs, outputs)
configs_param = get_parametrisation_configs(inputs_a, inputs_b, outputs, outputs)

# dict to store all polytopes
all_inequivalent_polytopes = {}
number_all_visited_polytopes = {}

# Generate vertices and it's permutations for the bell polytope
vertices = get_deterministic_behaviors(inputs_a, inputs_b, outputs)
#vertices_param = np.array(
#    [parametrise_behavior(p, configs, configs_param, inputs_a, inputs_b, outputs, outputs) for p in vertices])
permutations_vertices = np.loadtxt('../data/relabels_dets/{}{}{}{}.gz'.format(4, 3, 2, 2)).astype(int)

# Generate Bell polytope -> actually it does not matter here if we use vertices or parametrised vertices
bell_polytope = ParamPolytope(vertices, permutations_vertices)


def get_all_face_classes(polytope, level=0, max_level=1):
    """Function that returns all classes of faces for a polytope"""
    # storage of the polytope if an equivalent was not yet found
    if level not in all_inequivalent_polytopes.keys():
        all_inequivalent_polytopes[level] = [polytope]
        number_all_visited_polytopes[level] = 1
    else:
        number_all_visited_polytopes[level] += 1
        if not polytope.equiv_under_initial(all_inequivalent_polytopes[level]):
            all_inequivalent_polytopes[level].append(polytope)

    # if max recursion depth is reached, use Double Description to get all classes
    if level >= max_level:
        return polytope.get_all_classes()

    # get a single face
    face = polytope.get_single_face()
    # generate classes
    classes = [face]
    new_classes = [face]
    while new_classes:
        c = new_classes.pop()
        for r in get_all_face_classes(c, level + 1, max_level=max_level):
            res = c.rotate(r)
            if not res.equiv_under_parent(classes):
                classes.append(res)
                new_classes.append(res)
                if level == 0:
                    print('number of classes: ', len(classes))
    return classes


classes = get_all_face_classes(bell_polytope, max_level=1)
