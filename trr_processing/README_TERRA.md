1 ./run_op_simulations pdb_name water_model
	
	Specify the name of the pdb file to be used (no file handle).
	Water model to be used for simulation.

The concentrations are chosen with arbitrary number of cosolvents to be 
inserted into the simulation run.
A copy of "run_mpi_op.sh" needs to be present in the same working directory.
It will run:
	sbatch run_mpi+op.sh pdb_name water number_of_cosolvents

2. run_mpi_op.sh
	./combine2gro.sh central_chamber water_box combined_gro
	./make_pull_mdp.sh file_for_pulling gro_file lenght number_of_cosolvent
	sbatch run_simul_terra.sh lenght number_of_cosolvents
	
	Output simul_eqil.*


4. analyze_simulations.sh
	call_calc_op.sh
		graph_calc_op.py
			Outputs: calculated_osmotic_pressure.dat
					molality molarity idean osmotic_pressure calculated_osmotic_pressure		
				graph.pdf --> Full force vs. time trajectory
				graph_skip.pdf ---> call_calc_op.sh skips a certain number of frames

	Next it summarizes the measurements for all concentrations done into an op_(Number of cosolvents)nmg.dat

	Next it creates statistics.dat
	
	stat.py
		

After this use the script in bps_2015 to visualize and you are done!

