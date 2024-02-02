# Network Analysis for System Biology

## Object Scripts
### network.py
### substrate.py
### rate.py: 
### interaction.py:
#### Summary
    The interaction class acts as a data store for parameters and rules that govern how two substrates in the network are affected.

#### Attributes
    identifier: str
        - a key unique to an interaction between two substrates and can be used to refer to an interaction
        - required
    affected: Substrate
        - substrate object of the substrate whose value will change in this interaction between two substrates
        - required
    effector: Substrate
        - substrate object of the substrate whose value will produce the change in the affected substrate's value
        - required
    rate: Rate
        - rate object whose value acts as the proportionality constant between the affected and effector substrates
        - required
    effect: int (-1, 1, None)
        - if effect is -1, then the relationship between the substrate effector and affected is inverse meaning that as the effector substrates grows, the affected substrate falls; vice versa for effect of 1; None should only be used if the rate between the two substrates is the intrinsic phophorylation or dephosphoylation rates of the affected substrate and there is no feedback mechanism at play (no dissociation and hill coefficients)
        - optional in special cases (see logic in __init__ function)
    dissociation: Rate
        - rate object that encodes the 50% dissociation rate of two substrates; include if the interaction follows a feedback regulation mechanism; this package utilizes a goodwin oscillator mechanism to establish positive and negative feedback relationships; for more information see: https://pubmed.ncbi.nlm.nih.gov/32212037/
        - optional
    hill_coefficient: Rate
        - rate object that encodes the multiplicity of the two substrates (number of binding sites required); include if the interaction follows a feedback regulation mechanism; this package utilizes a goodwin oscillator mechanism to establish positive and negative feedback relationships; for more information see: https://pubmed.ncbi.nlm.nih.gov/32212037/
        - optional

#### Methods
    __init__():
            initialize the appropriate interaction objects and inform of any exceptions
    update_rate():
            adjust the rate of the interaction; required for finding optimal network parameters
    get_interaction_paramters():
            obtain all associated rates for the interaction of interest

For more information, check doc strings of each function individually.

## Processing Scripts
- **parse.py**: Contains rate, substrate, interaction parsing functions. These will take in the associated csvs and creates lists of rate, substrate, and interaction objects for the network object. For information about how to use these functions please look into file.

## Evaluation Scripts
- **score.py**: This script contains a function that will calculate the unexplained error of the model. This function takes in the predictions, ground truths, list of substrate names that encode the order of the prediction columns and optional objective function. The function will by default calculate mean square error but can be changed if pass through a objective function with inputs of a single prediction and associated truth and output of a cost.

## Experimental Scripts
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
### Model Scripts
### Experimental Scripts

## Requirements
- numpy
- pandas
- matplotlib
- scipy
- tqdm
- joblib
- geneticalgorithm2
