# HiDimFD

## Generation of Interpretable Residuals for Fault Diagnosis based on Projection Techniques: Leveraging Variable Redundancy

**Abel Alberto Cuadrado Vega, Ignacio Díaz Blanco, José María Enguita González, Diego García Pérez and Ana González Muñiz**

*Paper submitted for publication in IEEE Transactions on Control Systems Technology*

---

### Requirements
[SR3 MATLAB Toolbox ](https://github.com/UW-AMO/sr3-matlab)

### Usage
- Simulations in Section III:

Execute:

		eftsloop

- Experiments in Section IV:

First select the type of fault (sporadic or random) by commenting out the line not corresponding from the following two:

		fault_type = SPORADIC_FAULT;
		fault_type = RANDOM_FAULT;
Then execute:

		
		sporadic_and_random_faults