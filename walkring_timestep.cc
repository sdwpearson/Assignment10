// 
// walkring_timestep.cc
//
// Time stepping module for 1d random walk on a ring
//

#include "walkring_timestep.h"
#include <random>
#include <iostream>
#include <mpi.h>

// Perform a single time step for the random walkers
//
// parameters:
//
//  walkerpositions: the positions of a number of walkers (note that
//                   the number of walker Z can be found from
//                   Z=walkerpositions.size())
//
//  N:               the number of possible positions. All positions
//                   in the walkerpositions array should remain
//                   between 0 and N-1
//
//  prob:            the probability to jump to the left. Also the
//                   probability to just right.  (the probability to
//                   stay on the same spot is thus 1-2p.)
//
// output:
//
//  the content of the walkerpositions arrays should have changed to
//  reflect the random movement of all walker (i.e., they will each
//  have been given a chance to move on position to the left or two
//  the right).
//
void walkring_timestep(rarray<int,1>& walkerpositions, int N, double prob)
{
    // Create MPI environment
    int size, rank;
    int Z = walkerpositions.size();
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int send_count = Z/rank;
    int recv_count = send_count;
    int root = 0;
    rarray<int,1> scattered_walkers(recv_count);

    //Distribute the array of walkers
    MPI_Scatter(walkerpositions.data(), send_count, MPI_INT, scattered_walkers.data(), recv_count,  MPI_INT, root,  MPI_COMM_WORLD);

    int seed = (rank+1)*13;
    static std::mt19937 engine(seed);
    static std::uniform_real_distribution<> uniform;

    std::cout << scattered_walkers;

    // move all walkers
    for (int i = 0; i < recv_count; i++) {
        double r = uniform(engine); // draws a random number
        if (r < prob) {
            // move to the right, respecting periodic boundaries
            scattered_walkers[i]++;
            if (scattered_walkers[i] == N)
                scattered_walkers[i] = 0;
        } else if (r < 2*prob) {
            // move to the left, respecting periodic boundaries
            if (scattered_walkers[i] == 0)
                scattered_walkers[i] = N-1;
            else
                scattered_walkers[i]--;
        } else {
            // walkerposition remains unchanged
        }
    }

    if(rank == root)
        MPI_Gather(scattered_walkers.data(), recv_count, MPI_INT, walkerpositions.data(), send_count, MPI_INT, root, MPI_COMM_WORLD);
    else 
        MPI_Gather(scattered_walkers.data(), recv_count, MPI_INT, NULL, send_count, MPI_INT, root, MPI_COMM_WORLD);
}


