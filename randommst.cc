#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>

using namespace std;

// find index of minimum element given a set of indices to check from
// remove the visited node (since now part of MST)?
int deletemin(vector<int> &unvisited, vector<float> &dists) {
    int ind = 0;                        // initialize to the first element in the indices vector in case others don't exist
    float min = dists[unvisited[ind]];       // initialize minimum 
    for (int i = 1; i < unvisited.size(); i++) {
        float curr = dists[unvisited[i]]; 
        if (curr < min) {
            min = curr;
            ind = i;
        }
    }

    int return_val = unvisited[ind]; // store parent node (to be removed from "unvisited" next)
    
    // remove now visited node
    iter_swap(unvisited.begin() + ind, unvisited.end() - 1);
    unvisited.pop_back();
 
    return return_val;
}

// generate vector of random reals between 0 and 1 of specified input length
vector<float> randfun(int randlen) {
    vector<float> randomvec;
    for (int i = 0; i < randlen; i++) {
        randomvec.push_back((float) (rand()) / (float) (RAND_MAX));
    }
    return randomvec;
}

// return vector of euclidian distances for dimension 2
vector<float> euclid2d(vector<float> &xcoords, vector<float> &ycoords, vector<int> &unvisited, int currparent) {
    vector<float> out_dist;
    float x0 = xcoords[currparent];
    float y0 = ycoords[currparent];
    for(int i = 0; i < unvisited.size(); i++) {
        // euclidean distance is square root of sum of squared differences between the two vectors
        float euclid_dist = sqrt(pow((xcoords[unvisited[i]]-x0),2) + pow(ycoords[unvisited[i]]-y0, 2));
        out_dist.push_back(euclid_dist);
    }
    return out_dist;
}

// return vector of euclidian distances for dimension 3
vector<float> euclid3d(vector<float> &xcoords, vector<float> &ycoords, vector<float> &zcoords, vector<int> &unvisited, int currparent) {
   vector<float> out_dist;
   float x0 = xcoords[currparent];
   float y0 = ycoords[currparent];
   float z0 = zcoords[currparent];
   for (int i = 0; i < unvisited.size(); i++) {
       float euclid_dist = sqrt(pow((xcoords[unvisited[i]]-x0),2) + pow(ycoords[unvisited[i]]-y0, 2)+pow((zcoords[unvisited[i]]-z0),2));
       out_dist.push_back(euclid_dist);
   }
   return out_dist;
}

// return vector of euclidian distances for dimension 4
vector<float> euclid4d(vector<float> &xcoords, vector<float> &ycoords, vector<float> &zcoords, vector<float> &wcoords, vector<int> &unvisited, int currparent) {
   vector<float> out_dist;
   float x0 = xcoords[currparent];
   float y0 = ycoords[currparent];
   float z0 = zcoords[currparent];
   float w0 = wcoords[currparent];
   for(int i = 0; i < unvisited.size(); i++) {
       float euclid_dist = sqrt(pow((xcoords[unvisited[i]]-x0),2) + pow(ycoords[unvisited[i]]-y0, 2) + pow((zcoords[unvisited[i]]-z0),2) + pow((wcoords[unvisited[i]]-w0),2));
       out_dist.push_back(euclid_dist);
   }
   return out_dist;
}

int main(int argc,char* argv[])
{
    // command line parsing //
    
    int numpoints;
    int numtrials;
    int dimension;
 
    // argv[0] is "./randmst" //
    if(argc != 5) {
        printf("You passed too few or too many command line arguments passed (you passed %d arguments!)\n\n", argc);
        return 1;
    }
    
    else if (argc == 5) {
        // convert arguments to integers //
        numpoints = atoi(argv[2]);
        numtrials = atoi(argv[3]);
        dimension = atoi(argv[4]);
    }

    // trial runs go to outfile
    // open a file in append mode ?
    ofstream outfile;
    outfile.open("trials.txt", ios::app);

    int trial_count = 0;
    float total_distance = 0;
    float max_on_trials = 0;

    // begin clock
    clock_t timebegin = clock();
    float timesec;

    while (trial_count < numtrials) {
        trial_count++;

        vector<float> distance;         // track shortest distances
        vector<int> parent;             // parent to travel from to get minimum distance
        vector<int> unvisited_nodes;    // not in MST yet
        int currparent = -1;            

        // make vector of unvisited nodes, set all distances to infinity / parents to null
        for (int i = 0; i < numpoints; i++) {
            unvisited_nodes.push_back(i);   
            distance.push_back(10);          // represents infinity
            parent.push_back(-1);           // respresents NIL 
        }

        // starting node
        distance[0] = 0;        // at some starting node, dist is 0 and parent still NIL

        if (dimension == 1) {
            while(unvisited_nodes.size() > 0) {
                currparent = deletemin(unvisited_nodes, distance); // find min distance (want to add to MST)

                // generate edges from the node we're at now to all other nodes (aside from the ones already in MST)
                vector<float> newedges = randfun(unvisited_nodes.size());
                for(int i = 0; i < unvisited_nodes.size(); i++) {
                    int currchild = unvisited_nodes[i];              // gives relevant index
                    if(distance[currchild] > newedges[i]) {          // update if distance smaller  
                        distance[currchild] = newedges[i];
                        parent[currchild] = currparent;              // update parent 
                    }
                }
            }
        }

        if (dimension == 2) {
            vector<float> xcoords = randfun(numpoints);
            vector<float> ycoords = randfun(numpoints);

          while (unvisited_nodes.size() > 0) {
            currparent = deletemin(unvisited_nodes, distance);   // find min distance (which we add to MST)
             vector<float> newedges = euclid2d(xcoords, ycoords, unvisited_nodes, currparent);
                for (int i = 0; i < unvisited_nodes.size(); i++) {
                    int currchild = unvisited_nodes[i];              // gives relevant index
                    if (distance[currchild] > newedges[i]) {          // update if distance smaller  
                        distance[currchild] = newedges[i];
                        parent[currchild] = currparent;              // update parent
                    }
                }
            }
        }

        if (dimension == 3) {
            vector<float> xcoords = randfun(numpoints);
            vector<float> ycoords = randfun(numpoints);
            vector<float> zcoords = randfun(numpoints);

            while (unvisited_nodes.size() > 0) {
                currparent = deletemin(unvisited_nodes, distance);  // find min distance (which we add to MST)
                vector<float> newedges = euclid3d(xcoords, ycoords, zcoords, unvisited_nodes, currparent);
                for (int i = 0; i < unvisited_nodes.size(); i++) {
                    int currchild = unvisited_nodes[i];              // gives relevant index
                    if(distance[currchild] > newedges[i]) {          // update if distance smaller  
                        distance[currchild] = newedges[i];
                        parent[currchild] = currparent;              // update parent
                    }
                }
            }
        }

        if (dimension == 4) {
            vector<float> xcoords = randfun(numpoints);
            vector<float> ycoords = randfun(numpoints);
            vector<float> zcoords = randfun(numpoints);
            vector<float> wcoords = randfun(numpoints);

             while (unvisited_nodes.size() > 0) {
                currparent = deletemin(unvisited_nodes, distance); // find min distance (which we add to MST)

                vector<float> newedges = euclid4d(xcoords, ycoords, zcoords, wcoords, unvisited_nodes, currparent);
                for (int i = 0; i < unvisited_nodes.size(); i++) {
                    int currchild = unvisited_nodes[i];              // gives relevant index
                    if(distance[currchild] > newedges[i]) {          // update if distance smaller  
                        distance[currchild] = newedges[i];
                        parent[currchild] = currparent;              // update parent
                    }
                }
            }
        }

        // end clock
        clock_t timeduration = clock() - timebegin;

        timesec = (float) timeduration / CLOCKS_PER_SEC;

        // collecting total distance
        for (int i = 0; i < distance.size(); i++) {
            total_distance += distance[i];
        }

        // updating max element
        if (*max_element(distance.begin(), distance.end()) > max_on_trials) {
            max_on_trials = *max_element(distance.begin(), distance.end());
        }
    }

    outfile << "** Avg CPU Time (sec), Max Edge, Avg Total MST Weight, Numpoints, Numtrials, Dimension **" << endl;
    outfile << "Avg CPU Time (sec): " << (double) timesec / (double) numtrials << endl;
    outfile << "Max Edge: " << max_on_trials <<  endl;
    outfile << "Avg Total MST Weight: " << total_distance / (float) numtrials << endl;
    outfile << "n (Numpoints): " << numpoints << endl;
    outfile << "Number of Trials (Numtrials): " << numtrials << endl;
    outfile << "Dimension: " << dimension << endl;
    outfile << "////////////////////////////////////////////////////////////////////////////////" << endl;
    outfile.close();

    cout << "Avg Total MST Weight: " << total_distance / (float) numtrials 
         << ", n (Numpoints): " << numpoints 
         << ", Number of Trials (Numtrials): " << numtrials 
         << ", Dimension: " << dimension << endl;
    
    return 0;
}
