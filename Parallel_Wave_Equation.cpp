#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <string>
#include <filesystem>

#include <stdlib.h>
#include <mpi.h>

// This version of the code is designed to run on a local machine with user visibility of std::cout or std::cerr channels console printing purposes

int i_max = 100; // setting the global number of rows
int j_max = 100; // setting the global number of columns

int tag_num = 1; // setting a global tag_num

int id; // the rank of the individual process
int p; // the total number of processes

// setting up the pointers used to determine the domain, used later in the code, initialised as nullptr
int* i_local = nullptr;
int* i_local_start = nullptr;
int* upper_neighbour = nullptr;
int* lower_neighbour = nullptr;

// setting up the pointers to arrays used to store information comunicated about decomposted domains
int* recv_i_local = new int[1]; // the local number of rows in each decomposed domain
int* recv_i_local_start = new int[1]; // the global starting position of each decomposed domain
int* recv_upper_neighbour = new int[1]; // the id of the neighbouring process above
int* recv_lower_neighbour = new int[1]; // the id of the neighbouring process below

// used as a pointer to the global gid, populated later in master_sets_intial_conditions function
double* global_grid = nullptr;

// used in all_do_iteration_local
double* local_past_grid = nullptr;
double* local_current_grid = nullptr;
double* local_future_grid = nullptr;

// creating two arrays to receive data from upper and lower neighbours
double* upper_recv = nullptr;
double* lower_recv = nullptr;

// creating two arrays to send data to upper and lower neighbours
double* upper_send = nullptr;
double* lower_send = nullptr;

// setting the global spatial parameters of the simulation 
double y_max = 10.0; // the maximum extent of the y axis
double x_max = 10.0; // the maximum extent of the x axis
double dx = x_max / ((double)i_max - 1); // setting the step size in the x direction
double dy = y_max / ((double)j_max - 1); // setting the step size in the y direction
double c = 1; // the wave speed

// setting the global temporal parameters of the simulation
double t_max = 10.0; // the total length of time for which the simulation runs
double t = 0.0; // the time at which the simulation begins
double dt = 0.1 * std::min(dx, dy) / c; // setting the size of the time steps

// setting the global iteration counting parameters used to determine when to write to file
int iteration = 1; // setting the iteration counter, begins at 1 with the first iteration
int it_out = 100; // setting at which iteration increment the output is written to file

std::string directory_string; // setting a global string used in directory creation, broadcast with master process id 0

// for developer testing purposes only, user calls are not recommended, details have been deliberately ommitted from readme.md in the repo
void global_grid_to_file(std::string input_string)
{
    std::stringstream fname;
    std::fstream f1;
    fname << "./test/global_grid_" << input_string << "_" << id << ".dat";
    f1.open(fname.str().c_str(), std::fstream::out);
    for (int i = 0; i < i_max; i++)
    {
        for (int j = 0; j < j_max; j++)
        {
            f1 <<  global_grid[j + i * j_max] << "\t";
        }
        f1 << std::endl;
    }
    f1.close();
}

// for developer testing purposes only, user calls are not recommended, details have been deliberately ommitted from readme.md in the repo
void full_grid_to_file(int iteration, int step)
{
    std::stringstream fname;
    std::fstream f1;
    fname << "./out/test/iter"<< "_" << iteration << "_p_" << id << "_step_" << step << ".dat";
    f1.open(fname.str().c_str(), std::fstream::out);
    for (int i = 0; i < (*recv_i_local + 2); i++)
    {
        for (int j = 0; j < (j_max + 2); j++)
        {
            f1 <<  local_current_grid[j + i * (j_max + 2)] << "\t";
        }
        f1 << std::endl;
    }
    f1.close();
}

// a function to get the current time and convert to a string
const std::string get_current_time()
{
    time_t now = time(0); // record the system time as variable now
    struct tm time_struct; // create a struct for holding current time
    char time_buffer[80]; // creating a buffer of chars to be populated with the system time
    time_struct = *localtime(&now); // populating time_struct
    strftime(time_buffer, sizeof(time_buffer), "%d%m%y_%H%M%S", &time_struct); // creating a string in ddmmyy_HHMMSS format

    return time_buffer;
}

// a function to create the necessary directory to store the output
void all_create_directory(std::string input_directory)
{
    std::filesystem::create_directories("./output/" + input_directory + "/"); // creating a directory within the output folder to store images
}

// a function to broadcast the current time used in output file creation according to the master process
std::string master_broadcasts_time()
{
    int char_buffer_size = 14;
    char* directory_buffer = new char[char_buffer_size];

    if(id == 0) // continuity of current time is ensured by limiting get_current_time call to only process id 0
    {
        std::string current_time_string = get_current_time(); // gets curren time as a string
        std::copy(current_time_string.begin(), current_time_string.end(), directory_buffer); // convert string to char array
    }

    // broadcasts the current time to the directory buffer
    MPI_Bcast(directory_buffer, 14, MPI_CHAR, 0, MPI_COMM_WORLD);

    directory_buffer[char_buffer_size - 1] = '\0'; // enforce null termination of the string
    std::string directory_string(directory_buffer); // convert char array back to a string
    all_create_directory(directory_string); // create a directory named as directory_string

    delete[] directory_buffer; // freeing storage
    return directory_string;
}

// a function to write the inner grid of a domain to a file, ignoring ghost cells
void inner_grid_to_file(int iteration, std::string directory_string)
{
    std::stringstream f_ig_name; // creating the stringstream
    std::fstream f_ig; // creating the file stream
    std::string it_string = std::to_string(iteration); // converting the iteration number to a string, used in file writing
    std::filesystem::create_directories("./output/" + directory_string + "/" + it_string + "/"); // creating the necessary directory to store iteration images
    f_ig_name << "./output/" + directory_string + "/" + it_string + "/p"<< "_" << id << "_iter_" << iteration << ".dat"; // writing the .dat file
    f_ig.open(f_ig_name.str().c_str(), std::fstream::out); // opening the file stream
    // a loop to populate the newly created .dat file with the data contained in local_current_grid
    for (int i = 1; i < (*recv_i_local + 1); i++)
    {
        for (int j = 1; j < (j_max + 1); j++)
        {
            f_ig <<  local_current_grid[j + i * (j_max + 2)] << "\t"; // iterating using a 1D array
        }
        f_ig << std::endl;
    }
    f_ig.close(); // closing the file
}

// a function to write the time taken by each process to complete main to a .txt file
void time_to_file(double duration, std::string directory_string)
{
    std::stringstream fname; // creating the stringstream
    std::fstream f1; // creating the file stream
    std::string duration_string = std::to_string(duration); // converting the duration to a string, used in file writing
    namespace fs = std::filesystem; // creating filesystem namespace
    fs::create_directories("./timings/" + directory_string + "/"); // creating the necessary directory to store timings
    fname << "./timings/" + directory_string + "/p"<< "_" << id << "_timing.txt"; // writing the timing file
    f1.open(fname.str().c_str(), std::fstream::out); // openign the file stream
    int iteration_input = iteration - 1;
    f1 << "Time taken by process " << id << " on a " << i_max << "x" << j_max << " grid with " << iteration_input << " iterations: " << duration << " seconds"; // writing timing text to file
    f1.close(); // closing the file
}

// a function to set initial conditions on global_grid
void master_sets_intial_conditions(double user_r=1.0, double user_x=3.0, double user_y=3.0)
{
    // a check to ensure only the master process calls master_sets_intial_conditions
    if(id != 0)
    {
        std::cerr << "Error has occurred, only process 0 can call master_sets_intial_conditions" << std::endl;
        exit(1);
    }

	double r_splash = user_r;
	double x_splash = user_x;
	double y_splash = user_y;

	for (int i = 0; i < i_max; i++)
    {
		for (int j = 0; j < j_max; j++)
		{
			double x = dx * i;
			double y = dy * j;

			double dist = sqrt(pow(x - x_splash, 2.0) + pow(y - y_splash, 2.0));

			if (dist < r_splash)
			{
				double h = 5.0*(cos(dist / r_splash * M_PI) + 1.0);

				global_grid[j + i * j_max] = h;
			}
		}
    }
}

// a function to get the master process to decide the information pertaining to all local domains, called by master process only
void master_sets_domain()
{
    // a check to ensure that only the master process calls master_sets_domain
    if(id != 0)
    {
        std::cerr << "Error has occurred, only process 0 can call master_sets_domain" << std::endl;
        exit(1);
    }

    // creation of a 1D array with i_max * j_max elements
    global_grid = new double[i_max * j_max] {0};

    // population of parameters corresponding to the local domains, to be scattered by master process
    i_local = new int[p]; // to hold the size of the decomposed domain for each process
    i_local_start = new int[p]; // to hold the global start of the decomposed domain for each process
    upper_neighbour = new int[p]; // to hold the id of the neighbour above for each process
    lower_neighbour = new int[p]; // to hold the id of the neighbour below for each process

    // setting up integer variables to reduce repeat operations within the upcoming loops
    int i_local_holder = i_max / p; // a holding variable, to be modified later depending on the outcome of decomposition
    int i_local_terminal = 0; // used in the determination of the size of the local domain of the final process

    // a logical statement to determine the i_local value for the final process
    if((i_max % p) == 0) 
    {
        // if the number of processes is a factor of the total number of rows the size of each domain is i_max / p
        i_local_terminal = i_local_holder;
    }
    else
    {
        // if not then number of rows in the final process is the remaining rows after all other process have been allotted rows on an i_max / p basis
        i_local_terminal = i_max - (i_max / p) * (p - 1); 
    }

    // looping over each process and populating the p sized arrays with values needed for communication
    for(int i = 0; i < p; i++)
    {
        if(i == 0) // a special case for the master process
        {
            i_local[i] = i_local_holder; // the master process always has a local domain with i_local_holder number of rows
            i_local_start[i] = 0; // the master process always begins at global row index 0
            upper_neighbour[i] = p - 1; // the master process always has an upper neighbour of the final process
            lower_neighbour[i] = i + 1; // the master process always has a lower neighbour of process id 1
        }
        else
        {
            if(i == p - 1) // a special case for the final process
            {
                i_local[i] = i_local_terminal; // the final process always has a local domain with i_terminal number of rows
                i_local_start[i] = i * (i_max / p); // the final process always has a starting global index of the sum of all rows of all previous processes
                upper_neighbour[i] = i - 1; // the same logic as all non-master processes, the neighbour above is final process id - 1
                lower_neighbour[i] = 0; // the final process always has a lower neighbour of process id 0
            }
            else // for all processes that are neither the master or final process
            {
                i_local[i] = i_local_holder; // all other processes have a local domain with i_local_holder number of rows
                i_local_start[i] = i * i_local[i]; // all other processes have their starting global index at their id * their i_local value
                upper_neighbour[i] = i - 1; // all other processes have an upper neighbour as their id - 1
                lower_neighbour[i] = i + 1; // all other process have a lower neighbour as their id + 1
            }
        }
    }
}

// a function to scatter the information of the local domains, as determined by the master process, called by all processes
void all_scatter_domain()
{
    // MPI_Scatter to send i_local information
    MPI_Scatter(i_local, 1, MPI_INT, recv_i_local, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Scatter to send i_local_start information
    MPI_Scatter(i_local_start, 1, MPI_INT, recv_i_local_start, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Scatter to send upper_neighbour information
    MPI_Scatter(upper_neighbour, 1, MPI_INT, recv_upper_neighbour, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Scatter to send lower_neighbour information
    MPI_Scatter(lower_neighbour, 1, MPI_INT, recv_lower_neighbour, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

// a function to send the relevant portions of global_grid to each process, called by all processes
void all_divide_domain(double* global_grid)
{
    // setting up a request for all processes
    MPI_Request request_domain;

    local_current_grid = new double[(*recv_i_local + 2) * (j_max + 2)] {0}; // each process creates an array to store the values that will be sent by process 0
    
    // setting up a non-blocking receive to store their grid data once intial conditions have been set
    MPI_Irecv(local_current_grid, (*recv_i_local + 2) * (j_max + 2), MPI_DOUBLE, 0, tag_num, MPI_COMM_WORLD, &request_domain);

    if(id == 0) // this part of the function is performed only by process 0
    {
        for(int procs = 0; procs < p; procs++) // loop through all ids
        {
            // create an array of doubles to store the grid that will be sent
            double* send_local_domain = new double[(i_local[procs] + 2) * (j_max + 2)] {0};
            // values are initialised to 0 before population due to the stencil being cross shaped leaving corner elements untouched
     
            for(int i = 1; i < (i_local[procs] + 1); i++) // process id 0 uses i_local[procs] to determine the number of rows in the local domain of each process 
            {
                for(int j = 1; j < (j_max + 1); j++)
                {
                    {
                        // iterate through global_grid, using inforation about each processes local domain to determine what data from global_grid each process gets
                        send_local_domain[j + i * (j_max + 2)] = global_grid[(j - 1) + (i_local_start[procs] + i - 1) * j_max]; // process id 0 uses i_local_start[procs] to determine the start on global_grid of the local domain of each process
                    }
                }
            }

            // process 0 MPI_Ssend to ensure synchronisation in order to safely reuse send_local_domain buffer
            MPI_Ssend(send_local_domain, (i_local[procs] + 2) * (j_max + 2), MPI_DOUBLE, procs, tag_num, MPI_COMM_WORLD);
            delete[] send_local_domain; // freeing newly created storage
        }
    }

    // a call to MPI_Wait to ensure no race conditions occur
    MPI_Wait(&request_domain, MPI_STATUS_IGNORE);
    
    // initialisation of past and future grids used when all_do_iteration is called here to avoid replication in iteratively called all_do_iteration_local function
    local_past_grid = new double[(*recv_i_local + 2) * (j_max + 2)] {0};
    local_future_grid = new double[(*recv_i_local + 2) * (j_max + 2)] {0};

    // copying across values from local_current_grid to local_past_grid
    for(int i = 0; i < ((*recv_i_local + 2) * (j_max + 2)); i++)
    {
        local_past_grid[i] = local_current_grid[i];
    }
}

// a function to perform iterations on the local domain of each process
void all_do_iteration_local(int user_boundary_condition=0)
{
    // user_boundary_condition = 0 is dirichlet
    // user_boundary_condition = 1 is zero gradient
    // user_boundary_condition = 2 is periodic

    // creating two arrays to receive data from upper and lower neighbours
    upper_recv = new double[j_max] {0};
    lower_recv = new double[j_max] {0};

    // creating two arrays to send data to upper and lower neighbours
    upper_send = new double[j_max] {0};
    lower_send = new double[j_max] {0};

    // loop to populate the respective send arrays
    for(int i = 1; i < (*recv_i_local + 1); i++)
    {
        for(int j = 1; j < (j_max + 1); j++)
        {
            if(i == 1) // if on the first row of each domain, ignoring ghost cells, send this to the neighbour above
            {
                upper_send[j - 1] = local_current_grid[j + i * (j_max + 2)];
            }
            else // if on the last row of each domain, ignoring ghost cells, send this to the neighbour below
            {
                lower_send[j - 1] = local_current_grid[j + i * (j_max + 2)];
            }
        }
        if(i == 1) // forcing the loop to jump from row index 1 to row index -3, i is incremented again at the top of loop to reach index -2
        {
            i = *recv_i_local - 1;
        }
    }

    MPI_Datatype ghost_layer; // creating an MPI_Datatype for the ghost cells
    MPI_Type_contiguous(j_max, MPI_DOUBLE, &ghost_layer); // making use of the fact that data is stored in a 1D array, therefore contiguous data is guaranteed
    MPI_Type_commit(&ghost_layer); // committing the datatype

    // creating a request for upper and lower neighbours
    MPI_Request request_upper;
    MPI_Request request_lower;

    // each process performs two non-blocking sends
    MPI_Isend(upper_send, 1, ghost_layer, *recv_upper_neighbour, tag_num, MPI_COMM_WORLD, &request_upper);
    MPI_Isend(lower_send, 1, ghost_layer, *recv_lower_neighbour, tag_num, MPI_COMM_WORLD, &request_lower);

    // each process performs two non-blocking receives
    MPI_Irecv(upper_recv, 1, ghost_layer, *recv_upper_neighbour, tag_num, MPI_COMM_WORLD, &request_upper);
    MPI_Irecv(lower_recv, 1, ghost_layer, *recv_lower_neighbour, tag_num, MPI_COMM_WORLD, &request_lower);

    // using MPI_Wait to ensure no race conditions
    MPI_Wait(&request_upper, MPI_STATUS_IGNORE);
    MPI_Wait(&request_lower, MPI_STATUS_IGNORE);

	// copy the upper_recv and lower_recv values to the ghost rows of each domain
    for(int i = 0; i < (*recv_i_local + 2); i++)
    {
        for(int j = 1; j < (j_max + 1); j++)
        {
            if(i == 0)
            {
                if(id == 0) // a special case for process 0
                {
                    if(user_boundary_condition == 0) // if dirichlet
                    {
                        local_current_grid[j + i * (j_max + 2)] = 0;
                    }
                    if(user_boundary_condition == 1) // if zero-gradient
                    {
                        local_current_grid[j + i * (j_max + 2)] = local_current_grid[j + (i + 1) * (j_max + 2)];
                    }
                    if(user_boundary_condition == 2) // if periodic
                    {
                        local_current_grid[j + i * (j_max + 2)] = upper_recv[j - 1];
                    }
                }
                else // all other processes enter this loop and populate their local_current_grid with data received from process above
                {
                    local_current_grid[j + i * (j_max + 2)] = upper_recv[j - 1];
                }
            }
            else
            {
                if(id == (p - 1)) // a special case for the final process
                {
                    if(user_boundary_condition == 0) // if dirichlet
                    {
                        local_current_grid[j + i * (j_max + 2)] = 0;
                    }
                    if(user_boundary_condition == 1) // if zero-gradient
                    {
                        local_current_grid[j + i * (j_max + 2)] = local_current_grid[j + (i - 1) * (j_max + 2)];
                    }
                    if(user_boundary_condition == 2) // if periodic
                    {
                        local_current_grid[j + i * (j_max + 2)] = lower_recv[j - 1];
                    }
                }
                else // all other processes enter this loop and populate their local_current_grid with data received from process below
                {
                    local_current_grid[j + i * (j_max + 2)] = lower_recv[j - 1];
                }
            }
        }
        if(i == 0) // forcing the loop to jump from the first row to the last row
        {
            i = *recv_i_local;
        }
    }

    // freeing storage
    delete[] upper_recv;
    delete[] lower_recv;
    delete[] upper_send;
    delete[] lower_send;

    // a loop to populate the columns
    for(int j = 0; j < (j_max + 2); j++)
    {
        for(int i = 1; i < (*recv_i_local + 1); i++)
        {
            if(j == 0)
            {
                if(user_boundary_condition == 0) // if dirichlet
                {
                    local_current_grid[j + i * (j_max + 2)] = 0;
                }
                if(user_boundary_condition == 1) // if zero-gradient
                {
                    local_current_grid[j + i * (j_max + 2)] = local_current_grid[(j + 1) + i * (j_max + 2)]; 
                }
                if(user_boundary_condition == 2) // if periodic
                {
                    // filling out the left ghost cell column
                    local_current_grid[j + i * (j_max + 2)] = local_current_grid[j_max + i * (j_max + 2)];
                }
            }
            else
            {
                if(user_boundary_condition == 0) // if dirichlet
                {
                    local_current_grid[j + i * (j_max + 2)] = 0;
                }
                if(user_boundary_condition == 1) // if zero-gradient
                {
                    local_current_grid[j + i * (j_max + 2)] = local_current_grid[(j - 1) + i * (j_max + 2)]; 
                }
                if(user_boundary_condition == 2) // if periodic
                {
                    // filling out the right ghost cell column
                    local_current_grid[j + i * (j_max + 2)] = local_current_grid[1 + i  * (j_max + 2)];
                }
            }   
        }
        if(j == 0) // forcing the loop to jump from the first column to the second last, incremented again at the top of the loop
        {
            j = j_max;
        }
    }

    // now each process performs an iteration locally
    for(int i = 1; i < (*recv_i_local + 1); i++)
    {
        for(int j = 1; j < (j_max + 1); j++)
        {
            local_future_grid[j + i * (j_max + 2)] = pow(dt * c, 2.0) * ((local_current_grid[j + (i + 1) * (j_max + 2)] - 2.0 * local_current_grid[j + i * (j_max + 2)] + local_current_grid[j + (i - 1) * (j_max + 2)]) / pow(dx, 2.0) + (local_current_grid[(j + 1) + i * (j_max + 2)] - 2.0 * local_current_grid[j + i * (j_max + 2)] + local_current_grid[(j - 1) + i * (j_max + 2)]) / pow(dy, 2.0)) + 2.0 * local_current_grid[j + i * (j_max + 2)] - local_past_grid[j + i * (j_max + 2)]; 
        }
    }

    // freeing the storage originally assigned to local_past_grid as it is no longer needed
    delete[] local_past_grid;

    // creating temporary pointers to allow for the pointers to be swapped and cycle through grids
    auto temp_current_grid = local_current_grid;
    auto temp_future_grid = local_future_grid;

    local_past_grid = temp_current_grid; // past grid is now the former current grid
    local_current_grid = temp_future_grid; // current grid is now the former future grid
    local_future_grid = new double[(*recv_i_local + 2) * (j_max + 2)] {0}; // future grid is a new grid initialised to 0

    MPI_Type_free(&ghost_layer); // freeing the datatype after each iteration
}

// a function to free all storage
void all_free_storage()
{
    if(id == 0) // freeing storage created only by process 0
    {
        delete[] i_local;
        delete[] i_local_start;
        delete[] upper_neighbour;
        delete[] lower_neighbour;
        delete[] global_grid;
    }

    delete[] recv_i_local;
    delete[] recv_i_local_start;
    delete[] recv_upper_neighbour;
    delete[] recv_lower_neighbour;

    delete[] local_past_grid;
    delete[] local_current_grid;
    delete[] local_future_grid;
}

int main(int argc, char *argv[])
{
    // initialising MPI requirements
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

    double start = MPI_Wtime(); // recording the start time

    if(id == 0) // process id 0, set the domain and initial conditions
    {
        master_sets_domain();
        master_sets_intial_conditions();
    }

    // preventing race conditions so process id 0 can set domain and initial conditions before any communications are attempted
    MPI_Barrier(MPI_COMM_WORLD);

    all_scatter_domain();
    directory_string = master_broadcasts_time(); // all processes get the same string from process id 0 for directory creation
    all_divide_domain(global_grid); // all processes now divide global_grid

    // preventing race conditions to ensure domain has been adequately divided by all process before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    while (t < t_max) // a loop to iterate while t is less than t_max
	{
		all_do_iteration_local(2); // perform the iterations with user defined boundary conditions
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << "Process " << id << " iteration " << iteration << " complete " << "   time: " << t << std::endl; // print to console

        if((iteration % it_out) == 0) // write to file if iteration number satisifies writing criteria
        {
            inner_grid_to_file(iteration, directory_string);
        }

        t += dt; // increment t by dt
        iteration++; // increment the iteration counter
	}

    all_free_storage(); // freeing storage, this step is deliberately included in timings

    double end = MPI_Wtime(); // recording the end time
    double duration = end - start; // calculating the elapsed time
    std::cout << "Time taken by process " << id << ": " << duration << std::endl; // reading to console
    time_to_file(duration, directory_string); // writing time elapsed for each process to file

    // setting a barrier to combat a known bug on Macintosh machines related to overindexing of vader_segment
    // see https://github.com/open-mpi/ompi/issues/5798
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
}
