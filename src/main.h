//
//  main.h
//
//  Created by peter.fennell on 10/06/2014.
//  Copyright (c) 2014 peter.fennell. All rights reserved.
//

#ifndef main_h
#define main_h

// ALROGITHM FUNCTIONS
void initial_states(vector<int> &sigma, int X0, int N);
void increment_neighbours(int node, int k, vector<int> &neighbs, vector<int> &infected_neighbs);
void decrement_neighbours(int node, int k, vector<int> &neighbs, vector<int> &infected_neighbs);
void initialize_infected_neighbs(vector<int> &infected_neighbs, vector<vector<int> > &AdjList, vector<int> &degrees, vector<int> &sigma, int N);
void initialize_lists(vector<int> &infected_nodes, vector<vector<int> > &susceptible_nodes, vector<vector<int> > &susceptible_nodes_zero, vector<int> &positions, vector<int> &infected_neighbs, vector<int> &s_m, vector<int> &sigma, int N);
void remove_infected_node(int node, vector<int> &infected_nodes, vector<int> &positions);
void remove_susceptible_node(int node, int m, vector< vector<int> > &susceptible_nodes, vector<int> &positions);
void add_infected_node(int node, vector<int> &infected_nodes, vector<int> &positions);
void add_susceptible_node(int node, int m, vector< vector<int> > &susceptible_nodes, vector<int> &positions);
void decrement_susc_neighbs(int node, int k, vector<int> &neighbs, vector< vector<int> > &susceptible_nodes, vector<int> &infected_neighbs, vector<int> &s_m, vector<int> &positions, vector<int> &sigma);
void increment_susc_neighbs(int node, int k, vector<int> &neighbs, vector< vector<int> > &susceptible_nodes, vector<int> &infected_neighbs, vector<int> &s_m, vector<int> &positions, vector<int> &sigma);
double get_g_epsilon(double epsilon);

// MATRIX CONSTRUCTION FUNCTIONS
void GetAdjList_ER(vector<vector<int> > &AdjList, vector<int> &degrees, int N, double z, int &k_max);
void GetAdjList_CG(vector<vector<int> > &AdjList, vector<int> &degrees, int N);
int GetAdjList_custom(vector<vector<int> > &AdjList, vector<int> &degrees, int &N, const char* filename, int &k_max);


void initial_states(vector<int> &sigma, int X0, int N)
{
    fill(sigma.begin(), sigma.end(), 0);    // sets all elements of sigma to 0
    int counter = 0, node;
    
    while(counter < X0)
    {
        node = rand()%N;
        if(sigma[node] == 0)
        {
            counter ++;
            sigma[node] ++;
        }
    }
    
}


void increment_neighbours(int node, int k, vector<int> &neighbs, vector<int> &infected_neighbs)
{
    for(int i=0; i<k; i++)
    {
        infected_neighbs[neighbs[i]] ++;
    }
    
}


void decrement_neighbours(int node, int k, vector<int> &neighbs, vector<int> &infected_neighbs)
{
    for(int i=0; i<k; i++)
    {
        infected_neighbs[neighbs[i]] --;
    }
    
}

void initialize_infected_neighbs(vector<int> &infected_neighbs, vector<vector<int> > &AdjList, vector<int> &degrees, vector<int> &sigma, int N)
{
    fill(infected_neighbs.begin(), infected_neighbs.end(), 0);    // sets all elements of infected_neighbs to 0
    
    for(int i=0; i<N; i++)
    {
        // If node is in state +1, increment neighbours
        if(sigma[i]==1)
            increment_neighbours(i, degrees[i], AdjList[i], infected_neighbs);
    }
    
}

void initialize_lists(vector<int> &infected_nodes, vector<vector<int> > &susceptible_nodes, vector<vector<int> > &susceptible_nodes_zero, vector<int> &positions, vector<int> &infected_neighbs, vector<int> &s_m, vector<int> &sigma, int N)
{
    // Clear all the vectors
    fill(s_m.begin(), s_m.end(), 0);    // sets all elements of s_m to 0
    infected_nodes.clear();
    susceptible_nodes = susceptible_nodes_zero;
    
    for(int i=0; i<N; i++)
    {
        if(sigma[i]==1)
        {
            positions[i] = infected_nodes.size() ;
            infected_nodes.push_back(i);
        }
        else
        {
            positions[i] = susceptible_nodes[infected_neighbs[i]].size();
            susceptible_nodes[infected_neighbs[i]].push_back(i);
            s_m[infected_neighbs[i]] ++;
        }
    }

}

void remove_infected_node(int node, vector<int> &infected_nodes, vector<int> &positions)
{
    //Remove the node from the list by swapping the node with the last node in the list and removing the last node in the list. This will not affect the positioning of the other nodes once the position of the swapped last node is accounted for
    infected_nodes[positions[node]] = infected_nodes.back();    // Move last node in list to nodes position
    positions[infected_nodes.back()] = positions[node];          // account for this change in the position vector
    infected_nodes.pop_back();  // remove the last node from the list
}

void remove_susceptible_node(int node, int m, vector< vector<int> > &susceptible_nodes, vector<int> &positions)
{
    //Remove the node from the list by swapping the node with the last node in the list and removing the last node in the list. This will not affect the positioning of the other nodes once the position of the swapped last node is accounted for
    susceptible_nodes[m][positions[node]] = susceptible_nodes[m].back();    // Move last node in list to nodes position
    positions[susceptible_nodes[m].back()] = positions[node];          // account for this change in the position vector
    susceptible_nodes[m].pop_back();  // remove the last node from the list
}

void add_infected_node(int node, vector<int> &infected_nodes, vector<int> &positions)
{
    positions[node] = infected_nodes.size();
    infected_nodes.push_back(node);
    
}
void add_susceptible_node(int node, int m, vector< vector<int> > &susceptible_nodes, vector<int> &positions)
{
    positions[node] = susceptible_nodes[m].size();
    susceptible_nodes[m].push_back(node);
    
}

void decrement_susc_neighbs(int node, int k, vector<int> &neighbs, vector< vector<int> > &susceptible_nodes, vector<int> &infected_neighbs, vector<int> &s_m, vector<int> &positions, vector<int> &sigma)
{
    int neighb, m2;
    
    
    
    for(int j=0; j<k; j++)
    {
        neighb = neighbs[j];
     
        // Each neighbour of recovered node has one less infected neighbour
        infected_neighbs[neighb] --;
        
        // Change neighbourz from s_m class to s_{m-1} class
        if(sigma[neighb] == 0)
        {
            m2 = infected_neighbs[neighb];
            remove_susceptible_node(neighb, m2+1, susceptible_nodes, positions);
            s_m[m2+1] --;
            add_susceptible_node(neighb, m2, susceptible_nodes, positions);
            s_m[m2] ++;
        }
    }
}

void increment_susc_neighbs(int node, int k, vector<int> &neighbs, vector< vector<int> > &susceptible_nodes, vector<int> &infected_neighbs, vector<int> &s_m, vector<int> &positions, vector<int> &sigma)
{
    int neighb, m;
    
    for(int i=0; i<k; i++)
    {
        neighb = neighbs[i];
        
        // Each neighbour of infected node has one more infected neighbour
        infected_neighbs[neighb] ++;
        
        // Change neighbourz from s_m class to s_{m+1} class
        if(sigma[neighb] == 0)
        {
            m = infected_neighbs[neighb];
            remove_susceptible_node(neighb, m-1, susceptible_nodes, positions);
            s_m[m-1] --;
            add_susceptible_node(neighb, m, susceptible_nodes, positions);
            s_m[m] ++;
        }
    }
}

double get_g_epsilon(double epsilon)
{
    if(epsilon == pow(10.0, 1.0))
        return 11;
    if(epsilon == pow(10.0, 0.0))
        return 1.84181;
    if(epsilon == pow(10.0, -1.0))
        return 0.483183;
    if(epsilon == pow(10.0, -2.0))
        return 0.144835;
    if(epsilon == pow(10.0, -3.0))
        return 0.0450572;
    if(epsilon == pow(10.0, -4.0))
        return 0.0141755;
    if(epsilon == pow(10.0, -5.0))
        return 0.00447547;
    else
    {
        cout << "epsilon not in the form 10^n, where n is an integer" << endl;
        return -1;
    }
}

void GetAdjList_ER(vector<vector<int> > &AdjList, vector<int> &degrees, int N, double z, int &k_max)
{
    vector<vector<int> > adjlist(N);
    vector<int> degrees_temp(N,0);
    vector<int> neighbs;
    
    int m = round(z*N/2);
    int node1, node2, k;
    
    bool mult;
    
    // Assign edges untill there are none left
    while (m > 0)
    {
        // Pick two nodes
        node1 = rand()%N;
        node2 = rand()%N;
        if(node1 == node2)
            continue;
        
        
        k = degrees_temp[node1];
        neighbs = adjlist[node1];
        mult = false;
        for(int i=0; i<k; i++)
        {
            if(node2 == neighbs[i])
                mult = true;
        }
        
        if(mult == true)
            continue;
         

        
        else
        {
            adjlist[node1].push_back(node2);
            adjlist[node2].push_back(node1);
            degrees_temp[node1] ++;
            degrees_temp[node2] ++;
        }
        m--;
        
    }
    
    AdjList = adjlist;
    degrees = degrees_temp;
    
    // Find the Max degree
    k_max=0;
    for(int i=0; i<N; i++)
    {
        if(degrees[i] > k_max)
            k_max = degrees[i];
    }
    
    
}

void GetAdjList_CG(vector<vector<int> > &AdjList, vector<int> &degrees, int N)
{
    vector<vector<int> > adjlist(N);
    vector<int> adjlist_i, degrees_temp(N);
    
    // In a complete graph, every node is neighbour with every other node
    for(int i=0; i<N; i++)
    {
        for(int j1=0; j1<i; j1++)
            adjlist_i.push_back(j1);
        
        for(int j2=i+1; j2<N; j2++)
            adjlist_i.push_back(j2);
        
        degrees_temp[i]=N-1;
        
        adjlist[i] = adjlist_i;
        adjlist_i.clear();
        
    }

    AdjList = adjlist;
    degrees = degrees_temp;
    
    
}

int GetAdjList_custom(vector<vector<int> > &AdjList, vector<int> &degrees, int &N, const char* filename, int &k_max)
{
    
    // READ IN NETWORK
    // Network file is a list of edges of the form "node1 node2" i.e.
    //
    // 1 100
    // 1 1409
    // 2 3901
    // etc.
    vector<vector<int> > adjlist;
    vector<int> empty_vec;

    ifstream Afile;
    Afile.open(filename);
    
    
    // Alternative?
    
    
    if(!Afile.is_open())
    {
        cout << "Afile not found (please make sure that file_address is less than 100 characters)" << endl;
        return 0;
    }

    int temp, max_element;
    vector<int> list;   // read in all edges to this list and keep track of maximum element
    
    while(!Afile.eof())
    {
        Afile >> temp;
        if(temp > max_element)
            max_element = temp;
        list.push_back(temp);
    }
    
    if(list.size()%2 == 1)
    {
        cout << "ERROR: Odd number of nodes read in, make sure file is appropriately constructed and there is no white space at end of file" << endl;
        return 0;
    }
    
    
    vector<int> indices, indices_inverse(max_element+1,-1);   // indices is a function that give the node index of each node i=0,...,N-1 as they appear (possibly unordered from 0 to N-1) in the text file. index_inverse maps the node number from the text file to a node number between 0 and N-1 in the program
    
    
    int n = 0;  // count of number of nodes
    int node1, node2;
    vector<int>::iterator it=list.begin();
    while(it!=list.end())
    {
        node1 = *it;
        if (indices_inverse[node1]==-1) // If it is the first time the node has appeared
        {
            indices_inverse[node1] = n;
            indices.push_back(node1);
            adjlist.push_back(empty_vec);
            n ++;
        }
        
        ++it;
        node2 = *it;
        if (indices_inverse[node2]==-1)
        {
            indices_inverse[node2] = n;
            indices.push_back(node2);
            adjlist.push_back(empty_vec);
            n ++;
        }
        ++it;
        
        adjlist[indices_inverse[node1]].push_back(indices_inverse[node2]);
        adjlist[indices_inverse[node2]].push_back(indices_inverse[node1]);
    }
    
    vector<int> degrees_temp;
    int max_degree = 0, d;
    for(int i=0; i<n; i++)
    {
        d = (int) adjlist[i].size();
        degrees_temp.push_back(d);
        if(d > max_degree)
            max_degree = d;
    }
    
    AdjList = adjlist;
    degrees = degrees_temp;
    k_max = max_degree;
    N = n;
    
    return 1;
}


#endif
