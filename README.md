# PCMRPP v9.1

Authors: Shripad V Deshpande, Dr. Harikrishnan R. 

Symbiosis Institute of Technology, Pune Campus, Symbiosis International (Deemed University), Pune, India

Emails: deshpande.shripad2@gmail.com, dr.rhareish@gmail.com

POMDP file Creator for Mobile Robot Path Planning

This is a software package written in C for creating a pomdp file for mobile robot path planning. Partially observable Markov Decision Process (POMDP) is a well-proven mathematical framework for decision making in uncertain environment. The POMDP algorithm warks on a .pomdp file which lists down detailed states, observations, and state transition probabilities of the model. For accurately defining the states and observations in discrete space, the size of POMDP model becomes very big and creating such a file manually is a cumbersome task since it involves meticulously formulating the probabilities of huge number of state Transitions and Observations.
This software implements a novel algorithm to programmatically generate the POMDP model. THis is implemented in C on Linux. This version is limited to only two mobile robots moving in an open 2D space with possibility of collision between them. The software creates a pomdp file which can be passed through a pomdp solver (e.g. APPL Toolkit - https://bigbird.comp.nus.edu.sg/pmwiki/farm/appl/). The resultant policy file either in the form of graph or a set of alpha vectors are used by an executioner program to generate the robot path.

Configuration of the pomdp file is converted to pomdx format and The code is tested with examples of POMDP model with 1152 states with non-zero Probability spread over 10 observations.
Background Environment : Two mobile robots are assumed to be moving in a 2-D bounded field without any static obstacles. Each of the robot calls itself as R1 and the other robot as R2. Each robot fixes the frame of reference by assuming its own direction always to straight northward indicated by φ=900. The robot does all computations with respect to this frame of reference.  The state of the robot can be formulated using combination of parameters. [14] proposes that some of the parameters are fully deterministic whereas some are uncertain leading to the state with mixed observability. In this research work, the state of the robot is modelled using 5 parameters – its own speed (r1s), the angular direction of the other robot (α), the direction in which the other robot is moving (θ), the speed of the other robot (r2s) and the distance between the two robots (p).
