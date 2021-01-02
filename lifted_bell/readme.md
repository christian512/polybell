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

*  **joint_lifted_bell_extended_pr_box.py** - This implements the
    same idea as above. However here we use joint lifting, which 
   was an idea that we had during discussions about liftings.
   There was a counter example, that concluded that these joint liftings
   are not realizable or are preserving bell inequality properties.