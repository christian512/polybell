# Bell Inequalities from Linear Programming
In this repository, I try to explain the usage of linear programming, for finding Bell inequalities, for a given 
probability distribution. I follow the paper *Bell nonlocality* by *Brunner et. al. (2014)*, which can be found on the [*arXiv*](https://arxiv.org/abs/1303.2849).

The connection between Bell inequalities and linear programming interested me, so I though a short python demo would be
good to learn about the topic by myself. This repo is just what I've done to understand it better. Maybe it helps someone.

## Local behaviors and Bell inequalities

Bell inequalities can be used to distinguish between local and non-local probability distributions $p(ab|xy)$.
Here $x,y$ are the inputs to the measurement apparatus in a Bell scenario
and $a,b$ are the possible outcomes. The set of all possible probabilities $p = \{p(ab|xy)\}$ is called the **behavior**.
If there are $m$ different choices for the inputs $x,y$ and $\Delta$ unique outputs for each input, the behavior $p$ has
dimension $\Delta^2 m^2$.

Given such a behavior $p$, we want to identify if the correlations are local or not. Therefore we use, that the local
polytope $\mathcal{L}$ and the no-signalling polytope $\mathcal{NS}$ are both convex. We can then use the *hyperplane
seperation method* to see that a hyperplane that seperates a $p \notin \mathcal{L}$ from the local polytope $\mathcal{L}$.

If $\tilde{p} \notin \mathcal{L}$, then there exists an inequality

$$ s \cdot p = \sum_{abxy} s_{xy}^{ab} p(ab|xy) \leq S_k$$

that is fulfilled by all $p \in \mathcal{L}$, but not by $\tilde{p}$. This inequality is called the Bell inequality.

### Deterministic representation of local behaviors
In the case of a local behavior, a deterministic representation can be found. Therefore a local hidden variable $\lambda$
is used to assign an output to each input:

$$ \lambda = (a_1, ...,a_m; b_1, ..., b_m)$$

where $a_1$ is the output that is returned when the input $x=1$ is given. If there are $\Delta$ different outcomes,
we can find that there are $\Delta^m \cdot \Delta^m = \Delta^{2m}$ possible $\lambda$.

We can define $\Delta^{2m}$ different local deterministic behaviors $d_\lambda(ab|xy) \in \mathcal{L}$, which can
be used to describe **any** local behavior.

$$ d_\lambda(ab|xy) = 1 \quad \text{if } a=a_x \text{ and } b=b_y \quad\text{; otherwise } 0$$

The $d_\lambda$ form are unit vectors in the space $\mathbb{R}^{\Delta^{2m}}$ with only one non-zero entry. Hence any
local behavior $p$ can be written in the form
$$ p = \sum_\lambda q_\lambda d_\lambda \qquad \text{ where } \sum_\lambda q_\lambda = 1 \text{ and } q_\lambda \geq 0$$.

As all $d_\lambda$ are known, finding the weights $q_\lambda$ is a problem of linear programming.

### Linear problem
The problem of finding the weights can be written as

$$ \text{find } q_\lambda \text{ such that}$$ <br>
$$ \sum_\lambda q_\lambda d_\lambda = p \quad ; \quad \sum_\lambda q_\lambda = 1 \quad ; \quad q_\lambda \geq 0$$ <br>

However this would not be a optimization problem, but a satisfiability test. To transform it into an optimization problem,
we can introduce a **slack variable** $y$. The problem then reads as:

$$ \text{minimize } y \text{ such that }$$ <br>
$$ \sum_\lambda q_\lambda + y= 1 \quad ; \quad \sum_\lambda q_\lambda d_\lambda + y \cdot p= p$$ <br>
$$ q_\lambda \geq 0 \quad ; \quad y \geq 0 $$ <br>

Now this is the *primal form* of the problem and a trivial solution for $y = 1$ exists. Thus we can not tell if the behavior
$p$ was local or not, since the start configuration already satisfies the constraints. However we can look at the **dual
form** of the problem. This can be obtained by following some straightforward rules (found e.g. on Wikipedia). The dual
form then reads as:

$$\text{maximize } S = s\cdot p - S_l \text{ such that}$$ <br>
$$ s \cdot d_\lambda - S_l \leq 0 \quad \forall \lambda = 1,..., \Delta^{2m} $$ <br>
$$ s \cdot p - S_l \leq 1 $$

Here $s \in \mathbb{R}^{\Delta^{2m}}$ and $S_l \in \mathbb{R}$. If we now assume that $p$ is local, i.e. can be written
as convex combination of local deterministic behaviors $d_\lambda$, we find that $S \leq 0$. But we only required that
$S \leq 1$. Thus if $S > 0$ the behavior $p$ is non local. The **Bell inequality** is then:

$$ s \cdot p \leq S_l $$

Thus we can identify any local behavior by this solving this optimization problem. We can also directly say, if a behavior
is non-local.

### Implementation
We formulated the linear problem before and only need to solve it using linear programming. Such an optimizer is included
in the *SciPy* package. However we need to rewrite the problem, such that the *SciPy* solver can handle the problem.
The implementation of the problem can be found here: 