//
//  main.h
//
//  Created by peter.fennell on 10/06/2014.
//  Copyright (c) 2014 peter.fennell. All rights reserved.
//

#ifndef GA_h
#define GA_h

void GA_sim(vector< vector<int> > &AdjList, vector<int> &degrees, int k_max, 
int N, double beta, double mu, int X0, double tend, int n_intervals, double 
interval, int N_trials, vector<vector<int> > &X_vec) { 

    // VARIABLES
    double t, dt;           // time and time step
    vector<int> sigma(N);   // State Vector
    int X;                  // Number of Infected individuals
    int SI;                 // Number of links joining susceptible nodes to infected nodes
    vector<int> infected_neighbs(N);        // Number of infected neighbour of each node

    vector<int> infected_nodes;                         // list of the infected nodes
    vector<vector<int> > susceptible_nodes(k_max+1);    // list of the susceptible nodes with m infected neighbours
    vector<vector<int> > susceptible_nodes_zero(k_max+1);// list of the susceptible nodes with m infected neighbours uninitialized
    vector<int> positions(N);                           // position of each node in its correponding list
    vector<int> s_m(k_max+1,0);                         // Number of susceptible nodes with m infected neighbs
       
    // TEMPORARY VARIABLES
    int count;
    int node, m;
    double c, R;
    double s_m_sum;
    double mult1 = 1.0*(RAND_MAX-2)/pow((double)RAND_MAX,2.0);    // multiplicative constant used in generating uniform numbers from interval (0,1)
    double eps = 1.0/RAND_MAX;    // additive constant used in generating uniform numbers from the interval (0,1) (as opposed to [0,1]);
    vector<int> X_vec_ind(n_intervals + 1);

    
    for (int i=0; i<N_trials; i++) {
 
        // INITIALIZE SIMULATION
        // Ininitalize state vectors
        
        initial_states(sigma, X0, N);       // Radomly initialize state vector with X0 1s and rest 0s
        initialize_infected_neighbs(infected_neighbs, AdjList, degrees, sigma, N);
        initialize_lists(infected_nodes, susceptible_nodes, susceptible_nodes_zero, positions, infected_neighbs, s_m, sigma, N);
        
        X = X0;
        X_vec_ind[0] = X;
        SI = 0;
        for(int j=0; j<N; j++)                  // Calcuate inital number of SI links
            if(sigma[j] == 0)
                SI += infected_neighbs[j];

        
        
        // SIMULATION
        t = 0;
        count = 1;  // This is a variable used to record the intervals
        
        while(1)
        {
            
            // DRAW THE EXPONENTIALLY DISTRIBUTED TIME STEP
            c = rand()*mult1 + eps;     // c is a Uniform random number in the interval (0,1)
            R = X*mu + SI*beta;         // R is the total transition rate of the network
            dt = -log(c)/R;
            t += dt;
            
            
            // CHECK IF NEW TIME REQUIRES REDORDING OR IF FINISHED
            // If time has moves into one or more new intervals, then record time at them intervals
            while((t > count*interval)&&(count <= n_intervals))
            {
                X_vec_ind[count] = X;
                count ++;
            }
            
            if(t > tend)    // If time has elapsed
                break;
            

            
            // PICK A NODE, CHANGE ITS STATE AND UPDATE THE TRANSITION RATES
            c = 1.0*rand()/RAND_MAX;
            if(c < X*mu/R){
                // RECOVER AN INFECTED NODE
                
                
                // Pick an infected node at random
                node = infected_nodes[rand()%X];    // randomly choose an infected node
                m = infected_neighbs[node];         // number of infected neighbs of node
                sigma[node] = 0;                    // Change the state of the chosen node
                X--;                                // One less Infected node
                SI += 2*m - degrees[node];          // change in number of SI links as result of infected node changing
                
                
                // Remove node from infected list and add to susceptible list
                remove_infected_node(node, infected_nodes, positions);
                add_susceptible_node(node, m, susceptible_nodes, positions);

                // Increment number of susceptible nodes with m infected neighbours
                s_m[m] ++ ;

                // change the neighbours to other lists and adjust s_m
                decrement_susc_neighbs(node, degrees[node], AdjList[node], susceptible_nodes, infected_neighbs, s_m, positions, sigma);

            }
            else{
                // INFECT A SUSCEPTIBLE NODE
                
                
                // Pick an m class with probability proportional to the total number of SI links in that m class i.e. s_m[m]*m
                c = 1.0*rand()/RAND_MAX;        // Uniform random number in the interval [0,1]
                c *= SI;                        // Random number in the interval [0, N-X].
                s_m_sum = 0;
                for(int m1=0; m1<=k_max; m1++)
                {
                    s_m_sum += s_m[m1]*m1;
                    if(c <= s_m_sum)
                    {
                        // choose an node with m_infected neighbours at random
                        node = susceptible_nodes[m1][rand()%s_m[m1]];
                        break;
                    }
                }
                m = infected_neighbs[node];
                sigma[node] = 1;
                X ++;
                SI -= 2*m - degrees[node];   // change in number of SI links as result of susceptible node changing
                
                // Remove the node from susceptible list and add to infected list
                remove_susceptible_node(node, m, susceptible_nodes, positions);
                add_infected_node(node, infected_nodes, positions);
                
                // Decrement number of susceptible nodes with m infected neighbours
                s_m[m] -- ;
                
                // change the neighbours to other lists and adjust s_m
                increment_susc_neighbs(node, degrees[node], AdjList[node], susceptible_nodes, infected_neighbs, s_m, positions, sigma);
            }
            
            if(X == 0)      // If there are no further infected nodes
            {
                while (count <= n_intervals)
                {
                    X_vec_ind[count] = 0;
                    count ++;
                }
                break;
            }

            
            
        }//end of one time step
                
    }//end of one simulation
 
    X_vec.push_back(X_vec_ind);
}


void GA_sim_CG(vector< vector<int> > &AdjList, vector<int> &degrees, int 
k_max, int N, double beta, double mu, int X0, double tend, int n_intervals, 
double interval, int N_trials, vector<vector<int> > &X_vec) { 

    // VARIABLES
    double t, dt;           // time and time step
    int X;                  // Number of Infected individuals
       
    // TEMPORARY VARIABLES
    int count;
    double r, R;
    double mult1 = 1.0*(RAND_MAX-2)/pow((double)RAND_MAX,2.0);    // multiplicative constant used in generating uniform numbers from interval (0,1)
    double eps = 1.0/RAND_MAX;    // additive constant used in generating uniform numbers from the interval (0,1) (as opposed to [0,1]);
    vector<int> X_vec_ind(n_intervals+1);
    
    for (int i=0; i<N_trials; i++) 
    {
        
        t = 0;
        X = X0;
        
        X_vec_ind[0] = X;
        count = 1;
        
        while(1)
        {
 
            // DRAW TIME STEP
            r = rand()*mult1+eps;       // Random number between 0 and 1
            R = (beta*X)*(N-X)+mu*X;    // Total transition rate
            dt = -log(r)/R;             // random number from exponential distribution with parameter R
            
            t += dt;
            while((t > count*interval)&&(count <= n_intervals))
            {
                X_vec_ind[count] = X;
                count ++;
            }
            
            if(t > tend)
                break;
            
            
            // DRAW NODE
            r = 1.0*rand()/RAND_MAX;
            if( r <= X*mu/R )
            {
                // Infected node recovers
                X --;
            }
            else
            {
                // Susceptible node gets infected
                X ++;
            }
            
            
            if(X == 0)      // If there are no further infected nodes
            {
                while (count <= n_intervals)
                {
                    X_vec_ind[count] = 0;
                    count ++;
                }
                break;
            }
            
        }// end of one time step 
    
    X_vec.push_back(X_vec_ind);
                
    }// end of one simulation 
    
    
}


#endif
