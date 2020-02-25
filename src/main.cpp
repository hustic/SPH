/*
*   main.cpp
*   Purpose: call the classes for communication and plotting
*   to run and output game of life simulation using MPI
*/

#include "MpiComm.h"
#include "gnuplot.h"
#include <mpi.h>
#include <fstream>
#include <sstream>
#include <ctime>

using namespace std;

/**
 * game_of_life function takes the class object mp as input.
 * The class object contains the domain grid after decomposition.
 * the function should be called inside a iteration loop
 */
void game_of_life(MpiComm *mp, int &p_row, int &p_col);

/**
 * file_out funtion outputs the game of life grid of each processes into seperate files.
 */
void file_out(int &rank, domain *domains, int ite);

/**
 * plot_out function calls the gnuplot class which plots all processes
 * into one plot for each iteration.
 * 
 * The function directly output the plots into an gif animation.
 * 
 * The dele argument is set true by default, which deletes the
 * datafiles after creating the gif animation.
 */
void plot_out(MpiComm *mp, int &ite_max, bool dele);


// main to start the game of life simulation
int main(int argc, char *argv[])
{
    // initialise maximum number of iterations
    int ite_max = 201;

    // initialise dimension of 2D domain
    int row = 100;
    int col = 100;

    // initialise MPI parameters 
    int rank, p;

    // initialise start and end for timing the code with <ctime>
    double end, start;
    
    // initialise and assign MPI parameters
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    // A barrier before starting the timing to ensure accuracy
    MPI_Barrier(MPI_COMM_WORLD);
    // start timing
    start = MPI_Wtime();

    // starting mpi communication constructor
    // the constructor calls class methods from "domain.cpp" for domain decomposition
    // game of life grid values is randomly assigned inside the domain class method  
    auto *mp = new MpiComm(row, col);

    // the first communication method sends and receives the neighbours of the grid.
    // the neighbours are put into 4 vectors top[], bot[], left[] and right[]
    // the first communication also stores information about where the neighbours are located
    mp -> first_comm();

    // the rows and columns assigned to each processor after decomposition
    int p_row = mp->domains->p_row[rank];
    int p_col = mp->domains->p_col[rank];

    // output the initial domain into a file, iteration = 0
    file_out(rank, mp->domains , 0);

    // loop of game of life till maximum iterations
    for (int i = 1; i <= ite_max; i++)
    {
        // calls the game of life function each iteration
        game_of_life(mp, p_row, p_col);

        // output the updated domain into a new file
        file_out(rank, mp->domains , i);

        //barrier to make sure all the processes before communicating again 
        MPI_Barrier(MPI_COMM_WORLD);
        // communicate again using information obtained from first communication
        mp -> recomm();
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // end of iteration loop, use a barrier to make sure recorded time is accurate
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    
    // end MPI process
    delete mp;
    MPI_Finalize();

    // print out time spent to run the simulation
    if(rank == 0) cout << "time " << end - start << endl;
    
    // Plot the output data files into gif animation.
    // Decontructing the object will delete all the 
    // delete output data files if the final argument is true.
    if (rank == 0) plot_out(mp, ite_max, true);

    return 0;
}

void file_out(int &rank, domain *domains, int ite)
{
    // get coordinate system information
    int remain = domains->remain_p[0];
    int p_col = domains->p_col[rank];
    int p_row = domains->p_row[rank];
    int p_row_idx, p_col_idx;

    // use sstream and fstream library to output data files
	stringstream fname;
	fstream fp;
    // file names depend on rank and iteration number
	fname << "rank_" << rank << "_ite" << "_" << ite << ".dat";

    // open file
	fp.open(fname.str().c_str(), ios_base::out);

    // loop about the coordinate of each individual cell
    for (int i = domains->row_cord[rank] - p_row; i < domains->row_cord[rank]; i++)
    {
        // get the row cell index to get the value from individual process grids
        p_row_idx = i - (domains->row_cord[rank] - p_row);
        for (int j = domains->col_cord[rank] - p_col; j < domains->col_cord[rank]; j++)
        {
            // get the col cell index to get the value from individual process grids
            p_col_idx = j - (domains->col_cord[rank] - p_col);

            // output data into file, using the format of "row col value" and vertically.
            fp << double (i + 0.5) << " " << double (j + 0.5) << " " << domains->life_o[p_row_idx*p_col + p_col_idx] << "\n";
        }
        // leave a space after every new row for required formatting using gnuplot image plot
        fp << endl;
    }
    // close file
	fp.close();
}


void game_of_life(MpiComm *mp, int &p_row, int &p_col)
{
    // initialise parameters
    int temp;
    int rank = mp->rank;

    // assigning the top and bot neighbous into a new grid with two extra column and rows
    // top and bot neighbours contains the 4 corners of the grid, 2 in each 
    for(int i = 0; i < mp->top.size(); i++)
    {
        mp->domains->life_n[i] = mp->top[i];
        mp->domains->life_n[(p_row + 1)*(p_col + 2) + i] = mp->bot[i];
    }

    // assigning the left and right neighbous into a new grid with two extra column and rows
    for (int j = 0; j < mp->left.size(); j++)
    {
        mp->domains->life_n[(j + 1)*(p_col + 2)] = mp->left[j];
        mp->domains->life_n[(j + 2)*(p_col + 2) - 1] = mp->right[j];
    }

    // loop around row of each cell to determine live or die
    for (int j = 1; j < p_row + 1; j++)
    {
        // loop around col of each cell to determine live or die
        for (int k = 1; k < p_col + 1; k++)
        {
            temp = 0;
            // loop around row of neighbour cells
            for (int h = -1; h < 2; h++)
            {
                // loop around col of neighbour cells
                for (int g = -1; g < 2; g++)
                {
                    // the the number of cells alive
                    temp += mp->domains->life_n[(j + h)*(p_col + 2) + k + g];
                }
            }
            // only count the 8 neighbours surrounding the cell, taking away itsefl
            temp -= mp->domains->life_n[j*(p_col + 2) + k];

            // if there are bigger than 3 or less than 2 neighbours the cell must dead
            if (temp > 3 || temp < 2)
            {
                mp->domains->life_o[(j - 1)*p_col + k - 1] = 0;
            }
            // if there are 3 neighbours the cell must live
            else if (temp == 3)
            {
                mp->domains->life_o[(j - 1)*p_col + k - 1] = 1;
            }
        }
    }

    // copy the new life grid
    for (int i = 0; i < p_row; i++)
    {
        for (int j = 0; j < p_col; j ++)
        {
            mp->domains->life_n[(i + 1)*(p_col + 2) + j + 1] = mp->domains->life_o[i*p_col + j];
        }
    }
}

// function to call the gnuplot class for making the animation
void plot_out(MpiComm *mp, int &ite_max, bool dele = true)
{
    // only make animation on a sinlge processer
    if (mp->rank == 0)
    {
        auto *plt = new gnuplot(ite_max, mp->p, mp->row, mp->col, mp->rank, dele);
        plt -> plot();
        delete plt;
    }
}
