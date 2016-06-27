/* Ignore runnerup.sh */

run_trpcage.sh
	Starts the simulation
	
	The main thing is that it will run a set of simulations at room temperature, 
	it will then denature the protein and then quench it.

run.sh
	Checks the process of the simulation.
	
	sbatch terra_trpcage.sh		// Needs min.gro to proceed
	sbatch denature.sh		// Needs simul.gro to procedd
	sbatch production.sh		// Needs denature.gro to proceed


