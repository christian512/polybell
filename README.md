# Bell Inequalities from Linear Programming

*This is still work in progress*

In this repository, I try to explain the usage of linear programming, for finding Bell inequalities, for a given 
probability distribution. I follow the paper *Bell nonlocality* by *Brunner et. al. (2014)*, which can be found on the [*arXiv*](https://arxiv.org/abs/1303.2849).

The connection between Bell inequalities and linear programming interested me, so I though a short python demo would be
good to learn about the topic by myself. This repo is just what I've done to understand it better. Maybe it helps someone.

### Note for development
Just some notes I don't want to forget about. Maybe I have time for that at some point..

* Currently the local deterministic behaviors are created just from the list of inputs and outputs. However we have to 
define the local deterministic behaviors in the same way as the probability distribution. I.e. the entries at the same
indices must correspond to the same setting of inputs and outputs. Maybe you can find a nicer way to set this
