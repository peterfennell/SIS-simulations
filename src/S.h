//
//  main.h
//
//  Created by peter.fennell on 10/06/2014.
//  Copyright (c) 2014 peter.fennell. All rights reserved.
//

#ifndef S_h
#define S_h

int S_sim(vector< vector<int> > &AdjList, vector<int> &degrees, int k_max, int 
N, double beta, double mu, int X0, double tend, int n_intervals, double 
interval, int N_trials, vector<vector<int> > &X_vec) { 

    // if the maximum degree is very large then the maximum possible rate r_max = beta*k_max will result in a very small time step
    if(k_max > 10^2)
        cout << "WARNING:\nMaximum degree k_max is very large, will lead to inefficient simulations. Gillespie Algorithm simulations are preferable in this instance\n" << endl;

    // DYNAMICS
    double lambda_max;              // Maximum possible transition rate
    if(mu > k_max*beta)
        lambda_max = mu;
    else
        lambda_max = k_max*beta;
  
    
    // TIME (SET ACCORDING TO PAPER)
    // Time Step
    double epsilon = pow(10.0, -3.0);                  // epsilon as defined in section 3.B, has to be in form pow(10, n) where n is an integer between -5 and 1 (this integer range is just for the next function get_g_epsilon(epsilon) to work)
    double g_epsilon = get_g_epsilon(epsilon);         // Solution of Eqn (20) in section 3.B for given level of epsilon
    if(g_epsilon == -1)                     // Error (see get_g_epsilon() in Ssim.h)
        return 0;
    double dt = g_epsilon/lambda_max;
    // dt = 1;                                          // Uncomment to manually set dt
    cout << "Time step dt = " << dt << endl;
   
    
    // VARIABLES
    double t;               // time
    vector<int> sigma(N);   // State Vector
    int X;                  // Number of Infected and Susceptible individuals
    vector<int> infected_neighbs(N);
    vector<int> infected_neighbs_old(N);    // Number of Infected neighbours of each node at the start of a time step
    vector<int> X_vec_ind(n_intervals + 1);

 
    // TEMPORARY VARIABLES
    int count;
    double r;
  
    for (int i=0; i<N_trials; i++) {
        
        
        // INITIALIZE SIMULATION
        // Ininitalize state vectors
        
        initial_states(sigma, X0, N);    // Radomly initialize state vector with fraction rho_0 of 1s and rest 0s
        X = X0;
        X_vec_ind[0] = X;
        initialize_infected_neighbs(infected_neighbs, AdjList, degrees, sigma, N);
    
        // SIMULATION
        
        t = 0;
        count = 1;  // This is a variable used to record the intervals
        
        
        while(1)
        {
            
            t += dt;    // Increment time
            infected_neighbs_old = infected_neighbs;    // In S updating, the probability of a susceptible node updating depends on the number of infected neighbours at the start of the time step (infected_neighbs_old)

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
            
            
            
            
            for(int node = 0; node < N; node ++){
                
            
                // UPDATE NODE WITH CERTAIN PROBABLITY
                r = 1.0*rand()/RAND_MAX;
                if(sigma[node] == 1)
                {
                    if(r < mu*dt) // Infected node recovers
                    {
                        sigma[node] = 0;
                        X --;
                        decrement_neighbours(node, degrees[node], AdjList[node], infected_neighbs);  // Each of the neighbours of node now have one less infected neithgours
                    }
                }
                else
                {
                    if(r < infected_neighbs_old[node]*beta*dt) // Susceptible node gets Infected
                    {
                        sigma[node] = 1;
                        X ++;
                        increment_neighbours(node, degrees[node], AdjList[node], infected_neighbs);  // Each of the neighbours of node now have one less infected neithgours
                    }
                    
                    
                }
                
            }
                        
            
        //end of one time step
        }
                
        X_vec.push_back(X_vec_ind);
                
    //end of one simulation
    }
    
        
    return 1;

}


#endif

