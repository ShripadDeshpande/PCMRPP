# PCMRPP v9.1

Authors: Shripad V Deshpande, Dr. Harikrishnan R. 

Symbiosis Institute of Technology, Pune Campus, Symbiosis International (Deemed University), Pune, India

Emails: deshpande.shripad2@gmail.com, dr.rhareish@gmail.com

POMDP file Creator for Mobile Robot Path Planning

This is a software package written in C for creating a pomdp file for mobile robot path planning. Partially observable Markov Decision Process (POMDP) is a well-proven mathematical framework for decision making in uncertain environment. The POMDP algorithm warks on a .pomdp file which lists down detailed states, observations, and state transition probabilities of the model. For accurately defining the states and observations in discrete space, the size of POMDP model becomes very big and creating such a file manually is a cumbersome task since it involves meticulously formulating the probabilities of huge number of state Transitions and Observations.
This software implements a novel algorithm to programmatically generate the POMDP model. THis is implemented in C on Linux. This version is limited to only two mobile robots moving in an open 2D space with possibility of collision between them. The software creates a pomdp file which can be passed through a pomdp solver (e.g. APPL Toolkit - https://bigbird.comp.nus.edu.sg/pmwiki/farm/appl/). The resultant policy file either in the form of graph or a set of alpha vectors are used by an executioner program to generate the robot path.

Building executable from source :

            gcc -o pcmrpp pomdpV9-1.h utility.c pomdpCreateV9-1.c -lm

Running the software package :

            ./pcmrpp -v [0...4] -o [FILE]

**OPTIONS -**

-v : verbosity level of messages-

      0 Least verbose crisp messages - the execution will be faster
      4 Highest level of verbosity - the execution will be slower

-o : output pomdp filename. No assumption on the extension. The name will be taken as it is.

-h : Help about options

**Configurability**
Following parameters can be configured via the header file **pomdpV9-1.h**
Descritisation levels of angles alpha and theta (Refer the documentation for the context of these angles)

            #define N_DIV 8 //e.g. with this value, the angle will be descrtised to nearest multiple of 360/8 i.e. 45.
            
If the angles alpha and theta need different descritisation levels then following two lines will achieve that.

            #define ALPHA_MAX N_DIV
            #define THETA_MAX N_DIV
            
Following lines decide descritisation levels of other state parameters -

            #define R1SPEED_MAX 3
            #define R2SPEED_MAX 3
            #define PROXIMITY_MAX 2
            #define ACTION_MAX 6
            
The number of states and number of Observations would be same and equal to product of discrete levels of the 5 parameters-ALPHA_MAX, THETA_MAX, R!SPEED_MAX, R2SPEED_MAX, PROXIMITY_MAX. The size of Transition Matrix would be STATES_MAX * STATES_MAX * ACTIONS_MAX. The size of Observation Matrix would be STATES_MAX * OBSERVATIONS_MAX * ACTIONS_MAX.

The following line decides the span of probability distribution of Observations.

            #define MAX_OBS 10
            
e.g. MAX_OBS value of 10 indicates that the probability of the robot getting an observation on taking an action in a state is distributed across  ten observations. This controls the sparseness of the O matrix and in turn controls the number of alpha vectors generated by the pomdp solver. Higher the number of alpha vectors higher time of computation during execution because the belief state will have large number of entries and dot prodyuct of the belief satate vector with alpha vector will take time.
MAX_OBS value of 1 indicates state and observation is same which means it is MDP and no more partially observable.
