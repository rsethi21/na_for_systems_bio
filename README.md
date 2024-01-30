# Network Analysis for System Biology

## Summary

### Object Scripts

### Processing Scripts

### Evaluation Scripts

### Experimental Scripts
- **model_fitting_and_evaluation.py**:
	- This script creates a network object with the requested inputs and in return 
- **generate_samples.py**:
	- This script creates a network object with the requested inputs and in return creates a "n" number of samples with external stimuli and internal network states as features/attributes and phagocytosis/microsphere uptake rates as outputs.
	```
	python3 model_fitting_and_evaluation.py [-h] -s SUBSTRATES -r RATES -i INTERACTIONS [-f FITTING] [-a ARGUMENTS]
                                       [-o OUTPUT] [-p PARAMETERS] [-m MULTI] [-n NUMBER] [-t [TOPLOT ...]]	
	```
- **[your_requirements].py**:
	- 
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
