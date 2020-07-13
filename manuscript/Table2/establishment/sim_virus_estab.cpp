#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// g++ -std=c++11 sim_virus_time.cpp `pkg-config --libs gsl` command for compiling

using namespace std;

// Within-host parameters
#define R0 7.69             // R0 value
#define k 5.0           // eclipse phase: I_1 -> I_2
#define delta 0.595        // cell death
#define c 10.0           // virus clearance
#define p 11200             // continuous viral production
#define mu 0.001

// Initial variables
#define T0 40000        // initial number of target cells
#define V0 1            // initial number of infectious virions
#define E0 0           // initial number of eclipse phase cells
#define I0 0            // initial number of infected cells
#define W0 0			// initial number of non-infectious virions

double B = p/delta;         // burst mean
double beta = R0*c*delta/((mu*p-delta*R0)*(double)T0);    // virus infectivity

// Random number generation with Mersenne Twister
gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);

// Main code
int RUN(double,int,int);

int main(int argc,char *argv[])
{
    double eps = atof(argv[3]);                           // treatment efficacy
    int model = atoi(argv[2]);                  // 0 = burst model, 1 = continuous output model
    int scenario = atoi(argv[1]);               // 0 = burst size reduction, 1 = infectivity reduction
    int success = 0;                             // survival counter
    
    int r_ind=0;                                   // repetition index
    // int repeats = 100000;                           // number of repetitions
    
    while (success < 10000)
    {
        gsl_rng_set(r,r_ind);                    // setting the seed
        success += RUN(eps,scenario,model);         // updating survival counter (+1 if resistant strain survived)
        r_ind++;
        
        if (r_ind > 10000000)
        {
            break;
        }
    }
    
    return(0);
}

// Stochastic simulation
int RUN(double eps, int scenario, int model)
{
    int T, E, I, V, W;              // auxiliary variables (sensitive, resistant type)
    double t;                    // auxiliary variable (time) 
    int return_value = 0;        // return value (0 if V=0, 1 if V>1 at the end)
        
    // Initialization of the system
    t = 0.0;
    T = T0;
    E = E0;
    I = I0;
    V = V0;
    W = W0;
    
    // Simulation
    // Burst model
    if (model == 0)
    {
        while(( V>0 || E > 0 || I > 0 ) && V + W <= 2000)
        {
            // Update
            int update = 0;         // verification of update (while = 0 keep on searching for the index to update)
            int ind = 0;            // index for the transition to update
          
            // Transition rate vector
            double rates[5];
            
            // Reduction of burst size scenario
            if (scenario == 0)
            {
                rates[0] = beta*(double)T*(double)V;            // virus infecting cell
                rates[1] = k*(double)E;                         // leaving eclipse phase
                rates[2] = delta*(double)I;                     // cell death
                rates[3] = c*(double)V;                         // virus clearance
                rates[4] = c*(double)W;							// virus clearance
            }
                
            // Reduction of infectivity scenario
            if (scenario == 1)
            {
                rates[0] = beta*(1-eps)*(double)T*(double)V;    // virus infecting cell
                rates[1] = k*(double)E;                         // leaving eclipse phase
                rates[2] = delta*(double)I;                     // cell death
                rates[3] = c*(double)V;                         // virus clearance
                rates[4] = c*(double)W;							// virus clearance
            }   
            
            // Increase of virus clearance
            if (scenario == 2)
            {
                rates[0] = beta*(double)T*(double)V;    // virus infecting cell
                rates[1] = k*(double)E;                         // leaving eclipse phase
                rates[2] = delta*(double)I;                     // cell death
                rates[3] = c*(double)V/(1-eps);                         // virus clearance
                rates[4] = c*(double)W/(1-eps);					// virus clearance
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
                            V = V + gsl_ran_poisson(r,(1-eps)*mu*B);
                            W = W + gsl_ran_poisson(r,(1-eps)*(1-mu)*B);
                        }
                        
                        else
                        {
                            V = V + gsl_ran_poisson(r, mu*B);
                            W = W + gsl_ran_poisson(r, (1-mu)*B);
                        }
                    }                
                    
                    // virus death (infectious)
                    else if (ind == 3)
                    {
                        V = V-1;
                    }

                    // virus death (non-infectious)
                    else
                    {
                        W = W-1;
                    }
                }
               
                ind++;
            }
        }
    }
    
    // Continuous output model
    else
    {
        while(( V>0 || E > 0 || I > 0 ) && V + W <= 2000)
        {
            // Update
            int update = 0;         // verification of update (while = 0 keep on searching for the index to update)
            int ind = 0;            // index for the transition to update
          
            // Transition rate vector
            double rates[7];
            
            // Reduction of virus production scenario
            if (scenario == 0)
            {
                rates[0] = beta*(double)T*(double)V;            // virus infecting cell
                rates[1] = k*(double)E;                         // leaving eclipse phase
                rates[2] = delta*(double)I;                     // cell death
                rates[3] = mu*p*(1-eps)*(double)I;              // inf virus production
                rates[4] = c*(double)V;                         // virus clearance
                rates[5] = (1-mu)*p*(1-eps)*(double)I;			// non-inf virus production
                rates[6] = c*(double)W;							// virus clearance
            }
                
            // Reduction of infectivity scenario
            if (scenario == 1)
            {
                rates[0] = beta*(1-eps)*(double)T*(double)V;    // virus infecting cell
                rates[1] = k*(double)E;                         // leaving eclipse phase
                rates[2] = delta*(double)I;                     // cell death
                rates[3] = mu*p*(double)I;                         // virus production
                rates[4] = c*(double)V;                         // virus clearance
                rates[5] = (1-mu)*p*(double)I;
                rates[6] = c*(double)W;
            }     
            
            // Increase of virus clearance
            if (scenario == 2)
            {
                rates[0] = beta*(double)T*(double)V;    // virus infecting cell
                rates[1] = k*(double)E;                         // leaving eclipse phase
                rates[2] = delta*(double)I;                     // cell death
                rates[3] = mu*p*(double)I;                         // virus production
                rates[4] = c*(double)V/(1-eps);                 // virus clearance
                rates[5] = (1-mu)*p*(double)I;
                rates[6] = c*(double)W/(1-eps);
            }

	    if (scenario == 3)
	    {
			rates[0] = beta*(double)T*(double)V;
			rates[1] = k*(double)E;
			rates[2] = delta*(double)I/(1-eps);
			rates[3] = mu*p*(double)I;
			rates[4] = c*(double)V;
			rates[5] = (1-mu)*p*(double)I;
			rates[6] = c*(double)W;
	    }
                      
            // Draw two uniform random numbers for updating
            double rand1, dt;
            rand1 = gsl_ran_flat(r, 0.0, 1.0);  
            dt = gsl_ran_exponential(r, 1/accumulate(rates,rates+7,0.));  
            
            // Time update
            t += dt;
                
            // Population update
            while (update == 0)
            {
                if (rand1 < accumulate(rates,rates+ind+1,0.)/accumulate(rates,rates+7,0.))
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
                    else if (ind == 4)
                    {
                        V--;
                    }
                    
                    // non-inf virus production
                    else if (ind == 5)
                    {
                    	W++;	
                    }

                    // non-inf virus clearance
                    else
                    {
                        W--;
                    }
                }
               
                ind++;
            }      
        } 
    }
            
    if (V + W > 0)
    {
        return_value = 1;
        ofstream file ("esttime_LN_V0_1_eps_" + to_string(eps) + "_sc_" + to_string(scenario) + "_model_" + to_string(model) + ".txt", ios::app);   // file output, average
        file << t;  
        file << "\n";
        file.close();
    }
    
    return(return_value);
}
