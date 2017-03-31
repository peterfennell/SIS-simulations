//
//  main.cpp
//  SIS_AS
//
//  Created by peter.fennell on 15/04/2014.
//  Copyright (c) 2014 peter.fennell. All rights reserved.
//


#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <ctime>
#include <sstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

#include "main.h"
#include "GA.h"
#include "AS.h"
#include "S.h"


int main(int argc, const char * argv[])
{
    
    // ARGUMENTS PASSED TO FUNCTION

    int N, net, X0, N_trials, Xdat, EXdat, PXTdat, Sim;
    // N: Number of nodes; net: 1=ER, 2=Complete Graph; N_trials: number of simulations to realize; Xdat: 1 to save data on evolution of X for each simulation, 0 to not save; EXdat: 1 to save data on evolution of E[X] (average evolution of X over each simulation), 0 to not save; PXTdat: 1 to save data on evolution of P(x,t) for each simulation, 0 to not save;
    double mu, beta, tend, z;
    // mu and beta: SIS recovery and infections rates repsectively; tend: Simulation length; z: mean degree (only for ER network); X0: initial number of infected nodes;
    char filename[100];
    
    if(argc < 12)
    {
        cout << "Not enough arguments supplied, programme terminated" << endl;
        return 0;
    }
    if(argc > 13)
    {
        cout << "Too many arguments supplied, programme terminated" << endl;
        return 0;
    }
    sscanf(argv[1], "%d", &Sim);
    sscanf(argv[2], "%lf",  &beta);
    sscanf(argv[3], "%lf",  &mu);
    sscanf(argv[4], "%d",  &X0);
    sscanf(argv[5], "%d",  &net);
    if(net == 1)
    {
        sscanf(argv[6], "%d", &N);
        sscanf(argv[7], "%lf",  &z);
        sscanf(argv[8], "%lf",  &tend);
        sscanf(argv[9], "%d",  &N_trials);
        sscanf(argv[10], "%d",  &Xdat);
        sscanf(argv[11], "%d",  &EXdat);
        sscanf(argv[12], "%d",  &PXTdat);
    }
    if(net == 2)
    {
        sscanf(argv[6], "%d", &N);
        sscanf(argv[7], "%lf",  &tend);
        sscanf(argv[8], "%d",  &N_trials);
        sscanf(argv[9], "%d",  &Xdat);
        sscanf(argv[10], "%d",  &EXdat);
        sscanf(argv[11], "%d",  &PXTdat);
    }
    if(net == 3)
    {
        sscanf(argv[6], "%s",  filename);
        sscanf(argv[7], "%lf", &tend);
        sscanf(argv[8], "%d",  &N_trials);
        sscanf(argv[9], "%d",  &Xdat);
        sscanf(argv[10], "%d",  &EXdat);
        sscanf(argv[11], "%d",  &PXTdat);
    }
    if(!((net==1)||(net==2)||(net==3)))
    {
        cout << "Network paramter mis-specified, should be 1 (ER) or 2 (Complete Graph) or 3 (Custom)" << endl;
        return 0;
    }
    
  
    // NETWORK
    vector< vector< int > > AdjList;        // Adjacency List, AdjList[i] contains a list of the neighbours of node i
    vector< int > degrees;                  // degrees[i] is the number of neithbours of node i
    int k_max;                              // k_max is the maximum degree
    if(net == 1)
    {
        GetAdjList_ER(AdjList, degrees, N, z, k_max);
        cout << "Numerical simulation of SIS dynamics on an Erdos-Renyi Graphs with " << N << " nodes and mean degree " << z << endl;
    }
    if(net == 2)
    {
        GetAdjList_CG(AdjList, degrees, N);
        k_max = N-1;
        cout << "Numerical simulation of SIS dynamics on a Complete Graph with " << N << " nodes." << endl;
    }
    if(net == 3)
    {
        cout << "Reading in Edge List..." << endl;
        if(GetAdjList_custom(AdjList, degrees, N, filename, k_max) == 0)
            return 0;
        cout << "Numerical simulation of SIS dynamics on the graph from text file " << filename << " with " << N << " nodes." << endl;

    }
    
    // DYNAMICS
    cout << "Infection rate " << beta << ", recovery rate " << mu << ". Initial number " << X0 << " of infected nodes." << endl;
    if(X0 > N)
    {
        cout << "ERROR: X0 is too large (X0 > N)" << endl;
        return 1;
    }
   
    
    // TIME 
    // Sampling Intervals
    int n_intervals = 100;          // Number of times to record the values for output
    double interval = tend/(1.0*n_intervals);
    
    // OUTPUT VARIABLES
    vector<vector<int> > X_vec;
   
    srand((int) time(NULL));    // Random number generator seed
    
    // TIMING VARIABLES
    clock_t tick, tock;
    tick = clock();

    cout << "STARTING SIMULATION OF " << N_trials << " TRIALS, EACH OF LENGTH " << tend << "\n" << endl;
    
    if(Sim == 1)
    {
        if(net == 2)
            GA_sim_CG(AdjList, degrees, k_max, N, beta, mu, X0, tend, n_intervals, interval, N_trials, X_vec);
        else
            GA_sim(AdjList, degrees, k_max, N, beta, mu, X0, tend, n_intervals, interval, N_trials, X_vec);
    }
    if(Sim == 2)
    {
        AS_sim(AdjList, degrees, k_max, N, beta, mu, X0, tend, n_intervals, interval, N_trials, X_vec);
    }
    if(Sim == 3)
    {
        int S_sim_completion;
        S_sim_completion = S_sim(AdjList, degrees, k_max, N, beta, mu, X0, tend, n_intervals, interval, N_trials, X_vec);
        if(S_sim_completion == 0)
            return 0;
    }
    
    
     
    
    tock = clock();
    double sim_time = 1.0*(tock - tick)/CLOCKS_PER_SEC;
    
    cout << "SIMULATION COMPLETE\nSimulation time was " << sim_time  << " or " <<  sim_time/N_trials << " seconds per realisation" << endl;
    
    // SAVE OUTPUT
    
    string sim_name, name = "_Xdat_", graph, N_term = "_N", beta_term = "_b", mu_term = "_m" ,extension = ".txt";
    
    if(Sim == 1)
        sim_name = "data/GA";
    if(Sim == 2)
        sim_name = "data/AS";
    if(Sim == 3)
        sim_name = "data/S";
    
    if(net == 1)
        graph = "ER";
    if(net == 2)
        graph = "CG";
    if(net == 3)
    {
        string name(filename);
        graph = name;
    }
    
   
    // X
    if(Xdat == 1)
    {
        
        cout << "Saving X Data" << endl;
        
        ofstream outfile;
        stringstream sstm;
        sstm << sim_name << name << graph << N_term << N << mu_term << mu << beta_term << beta << extension;
        string filename = sstm.str();
        outfile.open(filename.c_str());
        if(!outfile.is_open())
        {
            cout << "ERROR: X outfile '" << filename << "' not opened. Please ensure that directory 'data' exists." << endl;
            return 1;
        }
        
        
        for(int i=0; i<=n_intervals; i++)
        {
            for(int j=0; j<N_trials; j++)
            {
                outfile << 1.0*X_vec[j][i]/N << "\t";
            }
            
            outfile << endl;
        }
        
    }
    
    // EX, the Expected value of X over time
    
    if(EXdat == 1)
    {
        
        cout << "Saving E[X] Data" << endl;
 
        vector<double> EX(n_intervals+1,0);     // Expected value of X at each interval time
        ofstream outfileEX;
        string nameEX = "_EXdat_";
        stringstream sstmEX;
        sstmEX << sim_name << nameEX << graph << N_term << N << mu_term << mu << beta_term << beta << extension;
        string filenameEX = sstmEX.str();
        outfileEX.open(filenameEX.c_str());
        if(!outfileEX.is_open())
        {
            cout << "ERROR: E[X] outfile '" << filenameEX << "' not opened. Please ensure that directory 'data' exists." << endl;
            return 1;
        }
        
        double EXt;    // Expected value of X at each time t
        for(int i=0; i<=n_intervals; i++)
        {
            EXt = 0;
            for(int j=0; j<N_trials; j++)
                EXt += X_vec[j][i];
            
            EXt /= 1.0*N_trials*N;
            outfileEX << EXt << "\n";
        }
        outfileEX << endl;
    }
    
    // p(x,t)
    
    if(PXTdat==1)
    {
        
        cout << "Saving P(x,t) Data" << endl;

        
        ofstream outfilePXT;
        string namePXT = "_PXTdat_";
        stringstream sstmPXT;
        sstmPXT << sim_name << namePXT << graph << N_term << N << mu_term << mu << beta_term << beta << extension;
        string filenamePXT = sstmPXT.str();
        outfilePXT.open(filenamePXT.c_str());
        if(!outfilePXT.is_open())
        {
            cout << "ERROR: P(x,t) outfile '" << filenamePXT << "' not opened. Please ensure that directory 'data' exists." << endl;
            return 1;
        }
        
        vector<double> pxt, pxt_0(N+1,0);  // pxt[j] is the probability that there are j infected nodes at time t
        
        for(int i=0; i<=n_intervals; i++)
        {
            pxt = pxt_0;
            for(int j1=0; j1<N_trials; j1++)
            {
                pxt[X_vec[j1][i]] += 1.0/N_trials;
            }
            
            for(int j2=0; j2<=N; j2++)
            {
                outfilePXT << pxt[j2] << "\t";
            }
            outfilePXT << endl;
            
        }
        
    }
    
    std::cout << "Done" << endl;
    
}
