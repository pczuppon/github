#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// g++ -std=c++11 sim_virus.cpp `pkg-config --libs gsl` command for compiling

using namespace std;

// Within-host parameters
#define R0 7.69             // R0 value
#define k 5.0           // eclipse phase: I_1 -> I_2
#define delta 0.595        // cell death
#define c 10.0           // virus clearance
#define p 112000             // continuous viral production
#define mu 0.001            // fraction of infectious virus produced

// Initial variables
#define T0 4000         // initial number of target cells (30ml resp tract * 2*10^5 cells /ml)
#define V0 1            // initial number of viruses
#define E0 0           // initial number of eclipse phase cells
#define I0 0            // initial number of infected cells

double B = mu*p/delta;         // burst mean
double beta = R0*c*delta/((mu*p-delta*R0)*(double)T0);    // virus infectivity

// Random number generation with Mersenne Twister
gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);

// Main code
int RUN(double,int,int);

int main(int argc,char *argv[])
{
    double eps = atof(argv[3]);                 // treatment efficacy
    int model = atoi(argv[2]);                  // 0 = burst model, 1 = continuous output model
    int scenario = atoi(argv[1]);               // 0 = burst size reduction, 1 = infectivity reduction, 2 = clearance increase
    int success = 0;                             // survival counter
    
    int r_ind;                                   // repetition index
    int repeats = 1;                           // number of repetitions
    
    for (r_ind=0; r_ind<repeats; r_ind++)
    {
        gsl_rng_set(r,r_ind);                    // setting the seed
        success += RUN(eps,scenario,model);         // updating survival counter (+1 if resistant strain survived)
    }
    
    double avg;                                  // empirical survival probability of R
    avg = (double)success/(double)repeats;
    
    ofstream file ("surv_HN_V0_1_sc_" + to_string(scenario) + "_model_" + to_string(model) + ".txt", ios::app);   // file output, average
    file << eps;
    file << ",";
    file << avg;  
    file << "\n";
    file.close();
    
    return(0);
}

// Stochastic simulation
int RUN(double eps, int scenario, int model)
{
    int T, E, I, V;              // auxiliary variables (sensitive, resistant type)
    double t;                    // auxiliary variable (time) 
    int return_value = 0;        // return value (0 if V=0, 1 if V>1 at the end)
        
    // Initialization of the system
    t = 0.0;
    T = T0;
    E = E0;
    I = I0;
    V = V0;
    
    // Simulation
    // Burst model
    if (model == 0)
    {
        while(( V>0 || E > 0 || I > 0 ) && V <= 500)
        {
            // Update
            int update = 0;         // verification of update (while = 0 keep on searching for the index to update)
            int ind = 0;            // index for the transition to update
          
            // Transition rate vector
            double rates[4];
            
            // Reduction of burst size scenario
            if (scenario == 0)
            {
                rates[0] = beta*(double)T*(double)V;            // virus infecting cell
                rates[1] = k*(double)E;                         // leaving eclipse phase
                rates[2] = delta*(double)I;                     // cell death
                rates[3] = c*(double)V;                         // virus clearance
            }
                
            // Reduction of infectivity scenario
            if (scenario == 1)
            {
                rates[0] = beta*(1-eps)*(double)T*(double)V;    // virus infecting cell
                rates[1] = k*(double)E;                         // leaving eclipse phase
                rates[2] = delta*(double)I;                     // cell death
                rates[3] = c*(double)V;                         // virus clearance
            }
            
            // Increase of virus clearance
            if (scenario == 2)
            {
                rates[0] = beta*(double)T*(double)V;    // virus infecting cell
                rates[1] = k*(double)E;                         // leaving eclipse phase
                rates[2] = delta*(double)I;                     // cell death
                rates[3] = c*(double)V/(1-eps);                 // virus clearance
            }     
                      
            // Draw two uniform random numbers for updating
            double rand1, dt;
            rand1 = gsl_ran_flat(r, 0.0, 1.0);  
            dt = gsl_ran_exponential(r, 1/accumulate(rates,rates+4,0.));  
            
            // Time update
            t += dt;
                
            // Population update
            while (update == 0)
            {
                if (rand1 < accumulate(rates,rates+ind+1,0.)/accumulate(rates,rates+4,0.))
                {
                    update = 1;
                    
                    // cell infection           
                    if (ind == 0)
                    {
                        V = V-1;
                        T = T-1;
                        E++;
                    }
                   
                    // eclipse phase
                    else if (ind == 1)
                    {
                        E = E-1;
                        I++;
                    }
                    
                    // cell death
                    else if (ind == 2)
                    {
                        I = I-1;
                        if (scenario == 0)
                        {
                            V = V + gsl_ran_poisson(r,(1-eps)*B);
                        }
                        
                        else
                        {
                            V = V + gsl_ran_poisson(r, B);
                        }
                    }                
                    
                    // virus death
                    else
                    {
                        V = V-1;
                    }
                    
                }
               
                ind++;
            }
        }
    }
    
    // Continuous output model
    else
    {
        while(( V>0 || E > 0 || I > 0 ) && V <= 500)
        {
            // Update
            int update = 0;         // verification of update (while = 0 keep on searching for the index to update)
            int ind = 0;            // index for the transition to update
          
            // Transition rate vector
            double rates[5];
            
            // Reduction of virus production scenario
            if (scenario == 0)
            {
                rates[0] = beta*(double)T*(double)V;            // virus infecting cell
                rates[1] = k*(double)E;                         // leaving eclipse phase
                rates[2] = delta*(double)I;                     // cell death
                rates[3] = mu*p*(1-eps)*(double)I;                 // virus production
                rates[4] = c*(double)V;                         // virus clearance
            }
                
            // Reduction of infectivity scenario
            if (scenario == 1)
            {
                rates[0] = beta*(1-eps)*(double)T*(double)V;    // virus infecting cell
                rates[1] = k*(double)E;                         // leaving eclipse phase
                rates[2] = delta*(double)I;                     // cell death
                rates[3] = mu*p*(double)I;                         // virus production
                rates[4] = c*(double)V;                         // virus clearance
            }     
            
            // Increase of virus clearance
            if (scenario == 2)
            {
                rates[0] = beta*(double)T*(double)V;    // virus infecting cell
                rates[1] = k*(double)E;                         // leaving eclipse phase
                rates[2] = delta*(double)I;                     // cell death
                rates[3] = mu*p*(double)I;                         // virus production
                rates[4] = c*(double)V/(1-eps);                 // virus clearance
            }  
                      
            // Draw two uniform random numbers for updating
            double rand1, dt;
            rand1 = gsl_ran_flat(r, 0.0, 1.0);  
            dt = gsl_ran_exponential(r, 1/accumulate(rates,rates+5,0.));  
            
            // Time update
            t += dt;
                
            // Population update
            while (update == 0)
            {
                if (rand1 < accumulate(rates,rates+ind+1,0.)/accumulate(rates,rates+5,0.))
                {
                    update = 1;
                    
                    // cell infection           
                    if (ind == 0)
                    {
                        V--;
                        T--;
                        E++;
                    }
                   
                    // eclipse phase
                    else if (ind == 1)
                    {
                        E--;
                        I++;
                    }
                    
                    // cell death
                    else if (ind == 2)
                    {
                        I--;
                    }           
                    
                    else if (ind == 3)
                    {
                        V++;
                    }     
                    
                    // virus death
                    else
                    {
                        V--;
                    }
                    
                }
               
                ind++;
            }      
        } 
    }
            
    if (V > 0)
    {
        return_value = 1;
    }
      
    
    return(return_value);
}
