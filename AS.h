//
//  main.h
//
//  Created by peter.fennell on 10/06/2014.
//  Copyright (c) 2014 peter.fennell. All rights reserved.
//

#ifndef AS_h
#define AS_h

void AS_sim(vector< vector<int> > &AdjList, vector<int> &degrees, int k_max, int 
N, double beta, double mu, int X0, double tend, int n_intervals, double 
interval, int N_trials, vector<vector<int> > &X_vec) { 


    // if the maximum degree is very large then the maximum possible rate r_max = beta*k_max will result in a very small time step
    if(k_max > 100)
        cout << "WARNING:\nMaximum degree k_max is very large, will lead to inefficient simulations. Gillespie Algorithm simulations are preferable in this instance\n" << endl;

    // DYNAMICS
    double lambda_max;              // Maximum possible transition rate
    if(mu > k_max*beta)
        lambda_max = mu;
    else
        lambda_max = k_max*beta;
  
    
    // TIME (SET ACCORDING TO PAPER)
    // Time Step
    double dt = 1.0/(N*lambda_max);
    cout << "Time step dt = " << dt << endl;
   
    // Transition probabilities
    double pmu = mu*dt*N;
    double pbeta = beta*dt*N;
    
    // VARIABLES
    double t;               // time
    vector<int> sigma(N);   // State Vector
    int X;                  // Number of Infected and Susceptible individuals
    vector<int> infected_neighbs(N);
    
    // TEMPORARY VARIABLES
    int count;
    int node;
    double r;
    vector<int> X_vec_ind(n_intervals + 1);
     
    for (int i=0; i<N_trials; i++) {
        
        
        // INITIALIZE SIMULATION
        // Ininitalize state vectors
        
        initial_states(sigma, X0, N);    // Radomly initialize state vector with fraction rho_0 of 1s and rest 0s
        
        X = X0;
        X_vec_ind[0] = X;
        initialize_infected_neighbs(infected_neighbs, AdjList, degrees, sigma, N);
        

        // Simulation
        
        t = 0;
        count = 1;  // This is a variable used to record the intervals
        
        
        while(1)
        {
            
            t += dt;    // Increment time
                    

            
            // If time has moves into one or more new intervals, then record time at them intervals
            while((t > count*interval)&&(count <= n_intervals))
            {
                X_vec_ind[count] = X;
                count ++;
            }
            
            if(t > tend)    // If time has elapsed
                break;
            
            if(X == 0)      // If there are no further infected nodes
            {
                while (count <= n_intervals)
                {
                    X_vec_ind[count] = 0;
                    count ++;
                }
                break;
            }
            

            // PICK NODE AT RANDOM
            node = rand()%N;
            
            // UPDATE NODE WITH CERTAIN PROBABLITY
            r = 1.0*rand()/RAND_MAX;
            if(sigma[node] == 1)
            {
                if(r < pmu) // Infected node recovers
                {
                    sigma[node] = 0;
                    X --;
                    decrement_neighbours(node, degrees[node], AdjList[node], infected_neighbs);  // Each of the neighbours of node now have one less infected neithgours
                }
            }
            else
            {
                if(r < 1 - pow(1-pbeta,infected_neighbs[node])) // Susceptible node gets Infected
                {
                    sigma[node] = 1;
                    X ++;
                    increment_neighbours(node, degrees[node], AdjList[node], infected_neighbs);  // Each of the neighbours of node now have one less infected neithgours
                }
                
            //end of one time step
            }
            
            

        
        
            
        //end of one simulation
        }        
               
        X_vec.push_back(X_vec_ind); 
        
    }
    

}


#endif
