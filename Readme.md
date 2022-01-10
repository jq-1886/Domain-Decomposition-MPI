<div align="center">
	<h1>ACSE-6 MPI Coursework</h1>
    	<h1>Parallelisation of the Wave Equation</h1>
	<h3>This software package demonstrates a parallel implementation of a simulation of the wave equation with horizontal domain decomposition on a global grid</h3>
</div>

***

# Contents
- [Instructions](#Instructions)
- [Languages](#Languages)
- [Parallel_Wave_Equation.cpp](#Parallel_Wave_Equation.cpp)
- [Post_Processing.py](#Post_Processing.py)
- [Images](#Images)
- [Timings](#Timings)
- [HPC](#HPC)

***

# Instructions

- To compile `Parallel_Wave_Equation.cpp` please use any recent gnu compiler, user may opt to implement optimsation flags.
- Note, master process refers to a process with id/rank of 0.

All simulations run by `Parallel_Wave_Equation.cpp`, with outputs in the `./Images/` and `./Timings/` directories have the following parameters:
- A 100x100 `global_grid`, (`i_max` and `j_max` of 100).
- A maximum extent in the x and y directions (`y_max` and `x_max`) of 10.
- A spatial step size (`dx` and `dy`) of `x_max / ((double)i_max - 1)` and `y_max / ((double)j_max - 1)`.
- A wavespeed (`c`) of 1 metre per second.
- A total simulation time (`t_max`) of 10 seconds.
- A time step (`dt`) of 0.1 * the minimum of either `dx` or `dy` all divided by `c`.
- A `it_out` value of 100, meaning output .dat files every 100 iterations.
- Performed with periodic boundary conditions.

All simulations run by `Parallel_Wave_Equation_HPC.cpp`, with outputs in the `./HPC/Images/` and `./HPC/Timings/` directories have the following parameters:
- A 1000x1000 `global_grid`, (`i_max` and `j_max` of 1000).
- A maximum extent in the x and y directions (`y_max` and `x_max`) of 10.
- A spatial step size (`dx` and `dy`) of `x_max / ((double)i_max - 1)` and `y_max / ((double)j_max - 1)`.
- A wavespeed (`c`) of 1 metre per second.
- A total simulation time (`t_max`) of 10 seconds.
- A time step (`dt`) of 0.1 * the minimum of either `dx` or `dy` all divided by `c`.
- A `it_out` value of 500, meaning output .dat files every 500 iterations.
- Performed with periodic boundary conditions.

***

# Languages

- Languages used are C++17 and Python 3.2.1.
- It is recommended that the user compiles or interprets .cpp or .py files with these version or later.
- Details of package requirements are in [Parallel_Wave_Equation.cpp](#Parallel_Wave_Equation.cpp) and [Post_Processing.py](#Post_Processing.py) sections.
- [Parallel_Wave_Equation.cpp](#Parallel_Wave_Equation.cpp) performs the parallel implementation of the wave equation simulation.
- [Post_Processing.py](#Post_Processing.py) performs the post processing on .dat file outputs of [Parallel_Wave_Equation.cpp](#Parallel_Wave_Equation.cpp).
- See relevant sections for more details.

***

# Parallel_Wave_Equation.cpp

This source file is used to perform a parallel implementation of a simulation of the wave equation. Please see `Coursework_Assignment.docx` for an explanation of the maths behind this code. The functions within the source code are detailed below.

#### const std::string get_current_time()
Function records current system time at the point at which it is called as a string and returns this string. The returned string is used in file creation to store output .dat and .txt files in the ./images/ and ./timings/ directories.

Function is not called directly in main, but is called indirectly when all processes call `master_broadcasts_time()`, where a logical statement results in only  master process calling it. 

#### void all_create_directory(std::string input_directory)
Function creates a directory ./output/`input_directory`/ in which simulation output .dat files are stored by each process.

Function is not called directly in main, but is called indirectly by all processes within `master_broadcasts_time()`.

#### std::string master_broadcasts_time()
Function results in master process broadcasting, through a call to `MPI_Bcast`, a char array that holds the system time at which master process called `get_current_time()`. The char array is then converted back to a string once the broadcast is complete. This string is then used by all processes to create a directory called ./output/`input_directory`/ where `input_directory` is the string broadcast by the call to `MPI_Bcast`, finally this string is returned.

All processes call this function directly in main.

#### void inner_grid_to_file(int iteration, std::string directory_string)
Function iterates through `local_current_grid` to write this array to a .dat file stored in ./output/`directory_string`/`it_string`/, where `directory_string` is the string returned by `master_broadcasts_times` and `it_string` is a string that holds the iteration number after conversion of `iteration` from type int to type string.

All processes call this function directly in main.

#### void time_to_file(double duration, std::string directory_string)
Function records time taken by each process to complete the code in main and writes to a .txt file. `duration` is generated by all process calls to `MPI_Wtime` in main and `directory_string` is the string returned by `master_broadcasts_time`. Function creates a directory ./output/`directory_string`/ where a .txt file is created with naming template `p_"0"_timing.txt`, where "0" is each processes id/rank. The time taken recorded in `duration` is then written to this file.

All processes call this function directly in main.

#### void master_sets_intial_conditions(double user_r=1.0, double user_x=3.0, double user_y=3.0)
Function sets the initial conditions on `global_grid` array before decomposition and communication. Initial conditions are set using a half-sinusoidal splash. All example runs of the code use the same initial conditions for comparative purposes.

Only master process directly calls this function in main.

#### void master_sets_domain()
Function performs initial stages of domain decomposition on `global_grid`, storing all dimensional information of all local domains in arrays which are subsequently scattered in an all process call to `all_scatter_domain`. Information generated is:
- `i_local` an array that holds the number of rows in each local domain.
- `i_local_start` an array that holds the starting index of each local domain on `global_grid`.
- `upper_neighbour` an array that holds the id of the neighbour above each process.
- `lower_neighbour` an array that holds the id of the neighbour below each process.

For all the above arrays, the index of each element refers to the id/rank of the process that element belongs to. For example `i_local[1]` holds the local number of rows in process 1's domain.

Only master process calls this function directly in main.

#### void all_scatter_domain()
Function calls `MPI_Scatter` to scatter all information pertaining to local domain dimensions contained in the arrays generated by `master_sets_domain`, where the ith element in each array corresponds to the data that belongs to the process with ith id/rank. The receive buffers for all processes in this funciton are the `recv_i_local`, `recv_i_local_start`, `recv_upper_neighbour` and `recv_lower_neighbour` arrays, which are all of size one and initialised to 0 as a global variables.

All processes call this function directly in main.

#### void all_divide_domain(double* global_grid)
Function performs non-blocking sends and receives so that master process can send all processes their respective local domains following the setting of initial conditions. This is done by master process entering an interative loop and performing the following:
- Creating a `send_local_domain` array for each process, with elements equal to the number of elements in each processes local domain.
- Population of the `send_local_domain` via element-wise copying from `global_grid`.
- A call to `MPI_Ssend` in order to implement synchronisation so as to preserve the `send_local_domain` buffer until safe reuse is possible.

Prior to master process entering the iterative loop, all processes create the `local_current_grid` array, used as the receive buffer in the round synchronised communications. `local_current_grid` contains (`i_max` + 2) * (`j_max` + 2) elements, so as to accomodate upper and lower rows of ghost cells.

The final step is to copy the values of `local_current_grid` to `local_past_grid` through element-wise iteration.

All processes call this function directly in main, though there is a block of code where a logical statement lets only master process enter.

#### void all_do_iteration_local(int user_boundary_condition=0)
Function makes each process perform iterations of the wave equation on their local domains. First 2 pairs of arrays are created by all processes via calls to `new`, `upper_recv` & `lower_recv` and `upper_send` & `lower_send`. The latter pair  being populated by each process with the first and last rows of their local domain (ignoring rows intended to store ghost cells), to be subsequently communicated to the neighbouring processes above and below, where the communicated data is used as those processes upper and lower rows of ghost cells. The former pair are used as receive buffers to store the information sent by the neighbouring processes above and below, which will contain data to be used as that processes upper and lower row of ghost cells.

Within this function a custom MPI datatype is used, `ghost_layer`, which is a 1D array that contains `j_max` elements and is stored using `MPI_Type_contiguous` to take advantage of the continguous memory form of 1D arrays. A call to `MPI_Wait` occurs after the round of non-blocking `MPI_Isend` and `MPI_Irecv` communications.

Once all necessary data have been received by each process, the data contained in `upper_recv` & `lower_recv` are copied to row index 0 and -1 of each processes `local_current_grid` to act as ghost cells.

Once this is completed, all processes free storage created to store `upper_recv` & `lower_recv` and `upper_send` & `lower_send`.

Columns are then populated according to the user specified bondary conditions where:
- a `user_boundary_condition` of 0 implements dirichlet boundary conditions at the edges of the global domain.
- a `user_boundary_condition` of 1 implements zero-gradient boundary conditions at the edges of the global domain.
- a `user_boundary_condition` of 2 implements periodic boundary conditions at the edges of the global domain.

The iterative calculations are then performed by an element-wise iterative traversal of `local_current_grid` and `local_past_grid` to populate `local_future_grid`. Pointers to these arrays are then cycled through to swap them so that:
- the next instance of `local_past_grid` has the memory address of the previous instance of `local_current_grid`.
- the next instance of `local_current_grid` has the memory address of the previous instance of `local_future_grid`.
- the next instance of `local_future_grid` is a newly created array of size (`i_max` + 2) * (`j_max` + 2) elements, initalised to 0, to be populated by the subsequent call to `all_do_iteration_local`.

The memory used to store the `local_past_grid` is immediately freed after completion of the iterative wave equation calculations, prior to the swapping of pointers so as to prevent memory leaks.

All processes call this function directly in main.

#### void all_free_storage()
Function frees all storage created by `Parallel_Wave_Equation.cpp` to prevent memory leaks.

All processes call this function directly in main, though there is a block of code where a logical statement lets only master process enter.

***

# Post_Processing.py

This Python script creates images and a .gif of the iterations output by `Parallel_Wave_Equation.cpp`. The script is desigend to be called once after each execution of `Parallel_Wave_Equation.cpp`.

The code works by first looking for the directory within `./output/` with the latest modification time. It then enters this subdirectory of `./output/` and then enters each directory within this subdirectory that holds all processes output .dat files. Within of these iteration subdirectories, a merged .dat file is created that contains in one file all processes output. Here, natural sorting is used by a call to the `natsorted` function of the `natsort` library so as to correctly used the applied naming convention of the output .dat files. The merged .dat file is then converted to a .png file and saved to an `./output/"latest_subdirectory"/images/` directory. Once this has been performed on all iterative steps, the .png files are then collated as a list, again using the `natsorted` function of the `natsort` library to create a .gif animation of the entire simulation performed by `Pararllel_Wave_Equation.cpp`.

***

# Images

This directory contains all images output by `Parallel_Wave_Equation.cpp`, example simulations have been run on the local machine used for development.

***

# Timings

This directory contains all timings output by `Parallel_Wave_Equation.cpp`, example simulations have been run on the local machine used for development.

***

# HPC

This directory contains the following directories:
- `./HPC/bin/` holds the compiled executable of `Parallel_Wave_Equation_HPC.cpp`, used on the Dug HPC to run the code, compiled with the command `mpicxx -std=gnu++17 -g Parallel_Wave_Equation_HPC.cpp -o Parallel_Wave_Equation_HPC`.
- `./HPC/logs/mpi_rjs/` holds the log files of each job.
- `./HPC/output/` holds the output .dat, .png and .gif files for each simulation run on the Dug HPC.
- `./HPC/timings/` holds the output. txt timings files for each simulation run on the Dug HPC.

`pwe.job` is the .job file used to submit all jobs to the Dug HPC.

Finally, within the `./HPC/` directory is `Parallel_Wave_Equation_HPC.cpp` which is the HPC version of `Parallel_Wave_Equation.cpp` with all writings to console removed and altered simulation parameters detailed in instructions.
