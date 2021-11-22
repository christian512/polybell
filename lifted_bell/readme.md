This directory contains codes for experimenting with
bell inequalities in a lifted form. Liftings are discussed
in many papers, so I don't explain them further here.

* **lifted_bell_extended_pr_box.py** - For a 3 outcome PR box
we check if there exists any bell inequality of lifted form
  that is violated by the PR box for different number of inputs. We can vary the efficiency 
  of the box, but we keep it slightly above threshold, which
  is the interesting setting here. 
  Result was that only for m = 2, 4 there are bell inequalities
  of lifted form. We checked for up to m = 10.

* **partial_lifted_bell_extended_pr_box.py** - Checks whether there exist 
partially lifted forms of bell inequalities for arbitrary input and 3 outputs.
  This is the case for 3 inputs, where both parties can have partial liftings.
  For the 6 inputs case only one party can have a partial lifting.
  This might indicate that for higher number of inputs there are no such symmetries.
  But we could check this on a machine with more computing power.
  
* **partial_lifted_bell_symmetric.py** - Checks whether there exist symmetric
bell inequalities for arbitrary number of inputs and 3 outputs. We did not find such Bell inequalities.
  Could check if the code is correct, but I don't think so.
  
* **analyse_bell_inequalities.ipynb** - Analyses the bell inequalities output by 
one of these scripts above.
  