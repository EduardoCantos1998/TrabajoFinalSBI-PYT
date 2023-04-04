# Theory
Protein–protein interactions (PPIs) are central to biological systems; and predicting the interacting residues is useful for constructing PPI networks, analysing mutations, drug design, drug discovery and improve annotation of protein funciton.
Experimental techniques commonly employed to determine the structure of protein complexes at atomic-scale resolution include X-ray crystallography nuclear mag- netic resonance (NMR) spectroscopy, and cryo-electron microscopy (cryo-EM). Information about interface residues can also be obtained by alanine scanning mutagenesis experiments or various footprinting experiments, such as hydrogen/deuterium exchange or hydroxy radical footprinting. As they may come useful they have the problem of still being currently expensive and low throughput. That's why in silico techniques could prove useful to fill the gap and determine 3D structure and interactions, specially given current available knowledge, larger datasets and that GPU acceleration have enabled the training of deeper neural network architectures. Broadly, there are 3 categories for PPI site prediction: a method focused in Machine Learning; Structure-based methods and Sequence-based methods.
 
![overview](overview.png?raw=true)

There are a number of actively developed machine learning frameworks. A popular choice for traditional ML is SciKit-Learn (https://scikit-learn.org/stable/). For deep learning frameworks, popular choices include: PyTorch (https://pytorch.org), TensorFlow (https://www.tensorflow.org), and Theano (http://deeplearning.net/software/theano/).

For our approach we decided to do Random Forest.

### Random Forest
Random Forests is an ensemble method that combines several individual classification trees in the following way: from the original sample several bootstrap samples are drawn, and an unpruned classification tree is fitted to each bootstrap sample. The feature selection for each split in the classification tree is conducted from a small random subset of predictor variables (features). From the complete forest the status of the response variable is predicted as an average or majority vote of the predictions of all trees.

### ML tool: scikit-learn
scikit-learn is a popular Python library used for machine learning tasks such as classification, regression, and clustering. It provides a range of tools and algorithms for data preprocessing, feature extraction, model selection, and model evaluation.


## Quality results
Regarding the efficiency of our tool we tried calculating the Root Mean Square Deviation (RMSD) between the predicted site and the actual site. RMSD is a measure of the difference between two sets of coordinates. In this case, we calculated the RMSD between the predicted site and the actual site using a software tool such as PyMOL or VMD. To do this, we aligned the predicted site with the actual site using a structural superposition algorithm. Then, calculated the RMSD between the aligned sets of coordinates. If the RMSD value was low (typically less than 2 Å), this indicated a good prediction. If the RMSD value was high, this indicated a poor prediction.  
In addition to RMSD, other metrics can also be used to evaluate ligand site predictions, such as the enrichment factor or the area under the receiver operating characteristic curve (AUC-ROC).

### DSSP
We used

## References