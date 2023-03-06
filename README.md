# Ariane-postprocessing
Python scripts to analyze and plot outputs from the particle tracking software Ariane.

[What is Ariane?](http://ariane.lagrangian.free.fr/ariane.html)

## Quantitative mode
About six million particles were initialised at the Strait of Gibraltar and followed back in time until they reached one of the control sections, in a maximum of 78 years. 

The control sections are defined in the **sections.txt** file as segments, forming vertical planes in the zonal or meridional directions, or horizontal planes, as in the case of the lid over the domain to trap particles entering the atmosphere. They are numbered sequentially, starting with section number 1, which corresponds to the initial section. Section number 2 corresponds to the Strait of Sicily, section number 3 to the Gulf of Lions and the number 4 to the Northern Tyrrhenian Sea. Ariane requires the specification of T-points indices for sections definition.

In order to better visualise the pathways of the water masses, three quantitative sub-experiments (one for each section) were carried out considering only the subset of particles that reached the section in question.The output files of the three sub-experiments are in the respective folders **/Lions**, **/Tyrrhenian**, **/Sicily**, while those concerning the complete experiment are in **/complete_experiment**.

In each folder there are two .nc files: **ariane_statistics_quantitative.nc** and **ariane_positions_quantitative.nc**. The former is used to calculate the horizontal streamfunction of the vertically integrated transport, while the latter file is used to calculate transit times and binned TS diagrams. In addition, each folder contains the **namelist** file for the parameters configuration required by Ariane for the experiment set-up, and the .txt file of statistics.  

The **mesh_mask_nemo.nc** file contains the information concerning the spatial domain of the reanalysis fields. 

For more detailed information on Ariane configuration, read the file **ariane_tutorial_sep08.pdf**.


