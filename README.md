ECS 129 Project 4; By: Munir Nur and Zaid Al Rakabi

Program to test how similar protein sequences are using the BL62 matrix and kernel method for comparing sequences developed by Smale et al. For further information on how the algorithm works and how we tested it read our paper: Comparing Protein Sequences using String Kernels. This is a work in progress.


(note: might need to run through ssh in CSIF computers or a computer that has Linux for accurate distance floating-point computation)

To compile:
	type "g++ kernel.cpp" into the command line


To run:
	type "./a.out bl62.txt 0.01 (sequence 1 filename) (sequence 2 filename)" into the command line
	
	(0.01 is our tested beta-value, for our unweighted String Kernel and bl62.txt is the BLOSUM scoring matrix)


