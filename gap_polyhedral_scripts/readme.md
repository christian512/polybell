# Files for Polyehdral 

### **generate_adm_polyhedral.py**
This python script generates a GAP code that can be run with the [**polyhedral package**](http://mathieudutour.altervista.org/Polyhedral/index.html) for GAP.
You can generate such a script for a Bell polytope by running:

```shell
python generate_adm_polyhedral.py 2 2 2 2
```

This creates a GAP script in the **scripts** directory. You can execute the resulting script by:

```shell
gap scripts/2222_polyhedral.g
```
Note that you need to install the **polyhedral package** before you can execute these scripts.
You can find install instructions [here](http://mathieudutour.altervista.org/Polyhedral/index.html).