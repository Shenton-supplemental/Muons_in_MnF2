# Supplemental material for _Covalent charge-neutral muon state in antiferromagnetic MnF<sub>2</sub>_

M. H. Dehn,<sup>1,2,3 </sup> R. Scheuermann,<sup>4</sup> J. K. Shenton,<sup>5,*</sup>  S. J. Blundell,<sup>6</sup> W. A. MacFarlane,<sup>2,3,7</sup> D. Prabhakaran,<sup>6</sup> A. Suter,<sup>4</sup> N. A. Spaldin<sup>5</sup>and R. F. Kiefl<sup>1,2,3,&sect;</sup>


<sup>1</sup>Department of Physics and Astronomy, University of British Columbia, Vancouver, BC V6T 1Z1, Canada    
<sup>2</sup>Stewart Blusson Quantum Matter Institute, University of British Columbia, Vancouver, BC V6T 1Z4, Canada    
<sup>3</sup><span style="font-variant:small-caps;">Triumf</span>, Vancouver, BC V6T 2A3, Canada   
<sup>4</sup>Laboratory for Muon Spectroscopy, Paul Scherrer Institute, Villigen AG, Switzerland    
<sup>5</sup>Department of Materials, ETH Zurich, CH-8093 Zürich, Switzerland   
<sup>6</sup>Oxford University Department of Physics, Clarendon Laboratory, Parks Road, Oxford OX1 3PU, United Kingdom   
<sup>7</sup>Department of Chemistry, University of British Columbia, Vancouver, BC, V6T 1Z1, Canada    
<sup>*</sup> For queries about the supplemental material in this repository contact [J. Kane Shenton](mailto:kane.shenton@stfc.ac.uk).


---
In these notebooks we provide supplemental material for our work on understanding the behaviour of positive muons in MnF<sub>2</sub>.


### Candidate sites
  Key <span style="font-variant:small-caps;">VASP</span> input and output files for each candidate configuration are available in the subdirectory `./candidate_configs`. A basic summary of the sites can be found in the Jupyter notebook: [Candidate_sites.ipynb](https://github.com/Shenton-supplemental/Muons_in_MnF2/blob/master/Candidate_sites.ipynb). There we also examine the displacement of the ions due to the muon, showing that our 3x3x4 supercell is converged. 


### Convergence
  We provide input and output files for our key tests of convergence with respect to plane-wave cutoff energy and k-point sampling density in the subdirectory: `./convergence`. These tests are summarised in the Jupyter notebook: [Convergence.ipynb](https://github.com/Shenton-supplemental/Muons_in_MnF2/blob/master/Convergence.ipynb).

### Hubbard U
  In the Jupyter notebook: [Hubbard_U.ipynb](https://github.com/Shenton-supplemental/Muons_in_MnF2/blob/master/Hubbard_U.ipynb) we investigate and summarise the effects of the choice of the Hubbard U correction and exchange-correlation functional on the crystal and electronic structure of MnF<sub>2</sub>. In addition, we investigate the effects of these on the computed <sup>19</sup>F hyperfine interactions. 


These Jupyter notebooks may be previewed on [github](https://github.com/Shenton-supplemental/Muons_in_MnF2) or via the [Jupyter notebook viewer](https://nbviewer.jupyter.org/github/Shenton-supplemental/Muons_in_MnF2). The latter sometimes does a better job of rendering the inline LaTeX and is therefore preferred.

Note that the vasprun.xml files are currently compressed to save space. These must be uncompressed before some of the notebooks will run. This can be done using the `gzip` command on UNIX-based systems. For example: `find . -name vasprun.xml.gz -exec gzip -d {} \;`
