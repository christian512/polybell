Scripts in this directory use PANDA to enumerate vertices,
and find Bell Inequalities for special cases.

* **test_panda_wo_pr.py** - Just test the PANDA code, which 
includes the equivalence check, that we also use in the implementation.
  
* **test_panda_pr.py** - Tries to find Bell Inequalities for 
PR-Boxes with 3 outputs. The script uses information about the 
  contributing deterministic points and thus should be the fastest 
  way that we currently have to find all such Bell Inequalities.