# HiDimFD

## Generation of Interpretable Residuals for Fault Diagnosis based on Projection Techniques: Leveraging Variable Redundancy

**Abel Alberto Cuadrado Vega, Ignacio Díaz Blanco, José María Enguita González, Diego García Pérez and Ana González Muñiz**

*Paper submitted for publication in IEEE Transactions on Control Systems Technology*

---

### Requirements
[SR3 MATLAB Toolbox ](https://github.com/UW-AMO/sr3-matlab)

### Dataset
[Dataicann: vibration and current data of an induction motor ](https://digibuo.uniovi.es/dspace/handle/10651/53461)

Download the **dataicann.zip** archive and place the contained **dataicann.mat** file in the directory with the scripts (only needed for `sporadic_and_random_faults`).

### Usage
- Simulations in Section III:

Execute:

		eftsloop

- Experiments in Section IV:

First select in the script `sporadic_and_random_faults` the type of fault (sporadic or random) by commenting out the line not corresponding from the following two, which can be found in said script:

		fault_type = SPORADIC_FAULT;
		fault_type = RANDOM_FAULT;
Then execute the script:

		
		sporadic_and_random_faults