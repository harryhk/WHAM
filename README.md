# WHAM
The project contains code for *Weighted Histogram Analysis Methods* for free energy calculations. This project has been down in year 2011 and if you want to use the code, please cite the following paper:   

**Kun Huang** and A.E. Garcia "Free energy of translocation an arginine-rich cell-penetrating peptide across a lipid bilayer suggests pore formation." Biophys. J. 104(2): 412-420. 

# Include
Three versions of WHAM.

	/root 
		/Wham_final1 : WHAM with both point methods and histogram method. 
		/Wham_final_openMp : WHAM with openMP support 
		/Wham_final_GPU: WHAM with GPU support 

# Dependencies:
A linux makefile is provide for each version of WHAM. To properly compile the program, you need:
* [GNU Scientific library (GSL)](http://www.gnu.org/software/gsl/).
* To compile WHAM with GPU support, you need to install cuda c compiler.


