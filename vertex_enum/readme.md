This file tells you about the content of the files and what they do.

## facet_enum_lrs.py
Doing some facet enumeration with lrs -> this is just for testing

## pr_box_local_weight.py
Calculating the local weight of the PR box and there by obtaining the deterministic behaviors
that contribute to the PR box. Then setting up a file for lrs which has equalities for the PR box and the contribution 
det. behaviors.

## vertex_cleaner.py
Runs through the vertices found by LRS and then checks for classes. So we can identify the number of classes with this 
script.

## vertex_enum.py
enumerates the vertices for a given number of inputs / outputs with only the deterministic points as constraints.

## vertex_enum.py
enumerates vertices for determnistic constraint and additional constraints using the PR box

## vertex_test_pr.py
Test the vertices found by LRS if the equalize the PR box condition PR * B == 1 and counts the number of vertices as well
as calculating the number of different classes that equalize these.

