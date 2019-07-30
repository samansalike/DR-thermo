# DR-thermo
The DR-thermo algorithm can be used to curate experimental thermodynamic data from different sources to generate reliable data-banks such that all measurements are modified to be thermodynamically consistent and any missing reaction or formation Gibbs energies can be predicted either directly or using group contributions. The algorithm can also be used to identify “unobservable” reaction and formation Gibbs energies, which can be targets for future experiments. MATLAB (preferably the latest version) is required to execute the code

RUNNING DR-Thermo 

To run DR-Thermo, save the reaction set in KEGG format and the relevant experimental values of reaction and formation energies into the "data" folder ( as well as pKa information and group decomposition matrix ) and run main.m from the "Main" folder. Refer to the Folder information given below for more information. 

FOLDER INFORMATION 

Data: 

The reaction set in KEGG format (representing the thermodynamic constraints) should be saved in reactions.txt file in the ‘reactions’ folder. Each reaction in this file is identified by the Reaction ID (RID) and each compound in the reaction is identified by the Compound ID (CID). The corresponding thermodynamic data for standard transformed reaction Gibbs free energies (dG’_r) and formation Gibbs free energies (dG’_f) corresponding to each RID and CID should be saved in the ‘reaction_dG.txt’ and ‘formation_dG.txt’ files respectively. When there are multiple Gibbs free energy measurements for the same reaction (or compound), either under the same or different reaction conditions, they should be included in the data with the same RID (or CID). The presence of such replicates would increase the redundancy of the DR algorithm. The pKa values should be saved as a pKa.mat file. The group decomposition matrix should be saved as G.mat.
 
Main:

Main.m is the file to be used to run DR-thermo. It generates all “chemical energy” measurements using the inverse Legendre transform and reads the reaction set to determine the constraint matrix (stoichiometric matrix S). From the available measurement data, it then determines the missing reaction and formation Gibbs energies. Thereafter, all available measurements are reconciled and the unavailable measurements are categorised as either observable or unobservable. For the unobservable variables, the group contribution estimates are generated using all available data and imputed onto the estimation step in the coaptation problem. Thereafter, all the reconciled data, in the form of standard Gibbs energies of reaction and formation corresponding to each reaction and compound in the RIDs and CIDs arrays, is compiled and stored in a “result” structure array. The RIDs and the CIDs corresponding the unobservable reactions and compound Gibbs energies are also displayed. When an unobservable variable cannot be estimated by group contributions, its resulting value in the result structure array will be represented by “NaN”. The RIDs and CIDs corresponding to these unobservable variables are stored in the RIDS_GC_NA and CIDS_GC_NA arrays. If these arrays are not present in the result, it means that group contribution estimates have been imputed for all unobservable variables. 

Misc:

The folder Misc contains different files that can be used to test the data reconciliation algorithm. Here, 2 reaction sets are given with the corresponding files for reactions.txt,  reaction_dG.txt, formation_dG.txt, G.mat, for those reaction sets. The files should be copied into the "data" folder to run the algorithm for the corresponding reaction set. The default files present in the "data" folder are from the 87 reaction set. 

87 reactions:

This reaction set of 87 reactions between 84 compounds is a subset of the reactions from the NIST-DECR database (as well as 13 additional redox potentials)  consisting of only those reactions which involve compounds whose formation energy values are available in the Alberty (2005) table of experimentally obtained Gibbs free energies of formation. This represents an example of a reaction set for which experimental data is completely available for all reaction and formation Gibbs energies (thus fully measured). 

473 reactions:

This reaction set of 473 reactions between 566 compounds consists of all reactions from the NIST-DECR database (as well as 13 additional redox potentials). The formation_dG.txt file here is obtained from the Alberty (2005) table and includes experimental formation energies for 225 compounds (117 of which are present in the 473 reactions). The G.mat file here is the group decomposition matrix for those 117 compounds. 


MATLAB FILE INFORMATION

Given below are discriptions for various ".m files" used in the repository. 


File: main.m

The main file to run data reconciliation for Gibbs free energy estimation. It reads the reactions from reaction_file (txt file with all the reactions in CIDS), measured apparent reaction and formation Gibbs free energies and generates the input data for recon_l.m. The apparent Gibbs energies are transformed into standard form using inverse Legendre transform ( see legendretransformF.m and legendretransfromR.m ) and the stoichiometric matrix is generated using parseKeggModel.m. Thereafter, the indices of the missing Gibbs energies are found by matching the available RIDS and CIDS from the measured Gibbs free energy data with the RIDS and CIDS from the reaction set (reaction_file). When there are replicates, the average of the values obtained after performing inverse Legendre transform is used. The variance of the values is used for constructing variance matrix. In the absence of replicates, the corresponding value in the variance matrix is unity. After performing reconciliation using recon_l.m, the results are sorted into the result.mat structure matrix which stores all the reconciled estimates, the unobservable CIDS and RIDS, and those not estimable using group contributions. 

Outputs:

	•	result.mat with vectors:
	⁃	dG0r_standard: The standard reaction Gibbs free energy estimates obtained from data reconciliation
	⁃	dG0f_Standard: The standard formation Gibbs free energy estimates obtained from data reconciliation
	⁃	CIDS_GC_NA: the CIDS corresponding to the compounds for which the modified group contribution method could not provide estimates ( due to insufficient measured data )
	⁃	RIDS_GC_NA: the RIDS corresponding to the reactions for which  the modified group contribution method could not provide estimates ( due to insufficient measured data )
	⁃	rids -  all RIDS present in the reaction set
	⁃	cids - all CIDS present in the reaction set

File: parseKeggModel.m (&  reaction2sparse.m)

Parses the reaction file with each line representing a reaction in KEGG format to obtain S, CIDS and RIDS

Inputs:

	•	ReactionStrings - cell array of each reaction in KEGG format
	•	arrow - arrow in the reaction (default is '=')

Outputs:

	•	S - Stoichiometric matrix
	•	cids - CIDS of all compounds in the reaction set
	•	rids - RIDS of all compounds in the reaction set

File: legendretransformF.m 

Applies the inverse Legendre transform to the transformed formation Gibbs energies. Maxstar.m and transform.m are used to perform the transform. 

Inputs:
 
	•	m - number of formation energies
	•	kegg_pKa - contains the relevant pKa, nH, charge etc information for each compound
	•	cids_m - CIDS of the formation energies that are to be inverse transformed
	•	pH - vector containing the pH values for each apparent Gibbs energy of formation measurement.
	•	T - vector containing the temperature values for each apparent Gibbs energy of formation measurement.
	•	I - vector containing the ionic strength values for each apparent Gibbs energy of formation measurement.
	•	dG0f_prime - the transformed Gibbs free energy of formation value that is to be converted to standard form ( or "chemical energy".

Outputs:

	•	dG0 - The measured standard Gibbs free energies of formation obtained after performing the inverse Legendre transform
	•	reverse_ddG0s - The difference between the standard  transformed ( or apparent) Gibbs free energies of formation and theGibbs free energies of formation.


File: legendretransformR.m 

Applies the inverse Legendre transform to the transformed reaction Gibbs energies. Maxstar.m and transform.m are used to perform the transform. 

Inputs:
	 
	•	S - stoichiometric matrix
	•	kegg_pKa - contains the relevant pKa, nH, charge etc information for each compound
	•	cids_m - CIDS of the formation energies that are to be inverse transformed
	•	pH - vector containing the pH values for each apparent Gibbs energy of reaction measurement.
	•	T - vector containing the temperature values for each apparent Gibbs energy of reaction measurement.
	•	I - vector containing the ionic strength values for each apparent Gibbs energy of reaction measurement.
	•	dG0r_prime - the transformed Gibbs free energy of reaction value that is to be converted to standard form 

Outputs:

	•	dG0 - The measured standard Gibbs free energies of reaction obtained after performing the inverse Legendre transform
	•	reverse_ddG0s - The difference between the standard  transformed ( or apparent) Gibbs free energies of reaction and theGibbs free energies of reaction.


File: recon_l.m

This file is used to perform the linear data reconciliation for Gibbs free energies. It calculates reconciled estimates for all measured Gibbs energies and computes prediction estimates for the unmeasured observable Gibbs energies. If prompted, it imputes the group contribution estimates for the unobservable Gibbs energies. It requires Gibbs free energy measurements, thermodynamic constraints (Stoichiometric matrix S), variance matrix (Sigma) and indices of the measured reaction and formation Gibbs energies corresponding to the columns and rows of S.  

Inputs: 
	•	 S - The (m x n) stoichiometric matrix for the n reactions between m compounds obtained from the reactions.txt file that represents the thermodynamic constraints for the available data 
	•	 y - vector of consisting of all measured reaction and formation Gibbs energies.
	•	index_rm -  indices of the "n" reaction Gibbs energies corresponding to the columns of S that are measured
	•	index_fm -  indices of the "m" formation Gibbs energies corresponding to the rows of S that are measured

Outputs:
	•	Recon_var= A vector of dimensions (m+n) x 1 representing the reconciled estimates. When the estimate is not available, the corresponding value is "NaN"
	•	index_unobservable= indices from the range 1:(m+n) that represent the unobservable Gibbs energies. 

File: Groupcontributions.m

This file is used to perform the modified group contribution method. By using all the available measured data, the thermodynamic constraints between them, estimates are obtained for all the group contributions using the constrained optimisation function fmincon. 

Inputs:

	•	A -  the reduced constraint matrix representing the thermodynamic constraints between the measured variables (obtained from the QR decomposition of S)
	•	G_r - measured reaction Gibbs energies
	•	G_f - measured formation Gibbs energies
	•	G - the group decomposition matrix
	•	index_fest - the indices of the unmeasured formation energies

File: costFunction.m
The optimisation function for fmincon in the group contribution method

Inputs:

	•	x - vector consisting of reconciled estimates of reaction, formation Gibbs energies and group contributions.
	•	b - vector consisting of measured reaction Gibbs energies
	•	Gf_m - vector consisting of measured formation Gibbs energies
	•	G_m - group decomposition matrix for the measured formation Gibbs energies

Outputs:

	•	J - The average sum squared error between measurements and estimates






