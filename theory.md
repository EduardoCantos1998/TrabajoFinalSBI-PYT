# Theory
Protein–protein interactions (PPIs) are central to biological systems; and predicting the interacting residues is useful for constructing PPI networks, analysing mutations, drug design, drug discovery and improve annotation of protein funciton.
Experimental techniques commonly employed to determine the structure of protein complexes at atomic-scale resolution include X-ray crystallography nuclear mag- netic resonance (NMR) spectroscopy, and cryo-electron microscopy (cryo-EM). Information about interface residues can also be obtained by alanine scanning mutagenesis experiments or various footprinting experiments, such as hydrogen/deuterium exchange or hydroxy radical footprinting. As they may come useful they have the problem of still being currently expensive and low throughput. That's why in silico techniques could prove useful to fill the gap and determine 3D structure and interactions, specially given current available knowledge, larger datasets and that GPU acceleration have enabled the training of deeper neural network architectures. Broadly, there are 3 categories for PPI site prediction: a method focused in Machine Learning; Structure-based methods and Sequence-based methods.
 
![overview](overview.png?raw=true)

For our approach we decided to do Random Forest. 

## Random Forest
Random Forests is an ensemble method that combines several individual classification trees in the following way: from the original sample several bootstrap samples are drawn, and an unpruned classification tree is fitted to each bootstrap sample. The feature selection for each split in the classification tree is conducted from a small random subset of predictor variables (features). From the complete forest the status of the response variable is predicted as an average or majority vote of the predictions of all trees.

## ML tool: scikit-learn
[scikit-learn](https://scikit-learn.org/stable/) is a popular Python library used for machine learning tasks such as classification, regression, and clustering. It provides a range of tools and algorithms for data preprocessing, feature extraction, model selection, and model evaluation. For deep learning frameworks, popular choices include: [PyTorch](https://pytorch.org), [TensorFlow](https://www.tensorflow.org), and [Theano](http://deeplearning.net/software/theano/).
## DSSP
To extract secundary structure information we used [DSSP](#references).


## [References](./README.md/#references)