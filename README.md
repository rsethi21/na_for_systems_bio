# Network Analysis for System Biology

## Summary

### Object Scripts

### Processing Scripts

### Evaluation Scripts

### Experimental Scripts
- **model_fitting_and_evaluation.py**:
	- This script creates a network object with the requested inputs and in return finds an optimal solution to network kinetic rates and compares the mean square error before and after fitting. This script also manually validates the model on other small datasets.
	- How to run:
        ```
        python3 model_fitting_and_evaluation.py [-h] -s SUBSTRATES -r RATES -i INTERACTIONS [-f FITTING] [-a ARGUMENTS]
                                       [-o OUTPUT] [-p PARAMETERS] [-m MULTI] [-n NUMBER] [-t [TOPLOT ...]]
        ```
	- To get information on these parameters:
        ```
        python3 model_fitting_and_evaluation.py -h
        ```
- **generate_samples.py**:
	- This script creates a network object with the requested inputs and in return creates a "n" number of samples with external stimuli and internal network states as features/attributes and phagocytosis/microsphere uptake rates as outputs.
	- How to run:
	```
	python3 generate_samples.py [-h] -s SUBSTRATES -r RATES -i INTERACTIONS [-o OUTPUT] [-p PARAMETERS] [-n NUMBER]
                           [-m MULTIPROCESS]	
	```
	- To get information on these parameters:
	```
	python3 generate_samples.py -h
	```
- **[your_requirements].py**:
	- You can design your own experiments. The network object has various toolings available to customize an analysis and generate data that can be graphed and manipulated for custom analysis.
## Inputs

## Outputs

## Requirements
- numpy
- pandas
- matplotlib
- scipy
- tqdm
- joblib
- geneticalgorithm2
