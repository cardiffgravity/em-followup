# em-followup
EM followup of Gravitational Wave sources

# Help

### Module Requirements in Python
Astropy 3.0.3
Requests 2.14.2

### Request_LCO
- This module contains a class for requesting observations with the [Las Cumbres Observatories](https://lco.global) network of telescopes.
- The proposal used in this porject only has access to the 0.4m telescopes, therefore this code is specific to those telescopes.

**Example of use:**

### download_data_LCO


**Example of use:**

### AAVSO_get_variables
- The Zooniverse project asks users to circle the variations observed in images taken 5-7 days apart.
- This module searches the AAVSO database using given RA and Dec limits and returns a list of coordinates for the known variable stars in 
that field.
- Circle are drawn around these known variables so that only "new" variables will be highlighted in the Zooniverse project.

**Example of use:**
```python
>>>fromRA = 9.1
>>>toRA = 10.8
>>>fromDec = -4.1 	
>>>toDec = 3.2

>>>RA, Dec = get_variables(fromRA, toRA, fromDec, toDec)

>>>print(RA)

9.11692
9.22075
9.32167
9.39625
9.45395

>>>print(Dec)

-0.52356
3.11739
1.96000
0.29417
-0.05333
```

### imcreate


**Example of use:**

### run_code

