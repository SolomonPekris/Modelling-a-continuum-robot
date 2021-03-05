# Modelling-a-continuum-robot
A dynamic model that calculates the position of a tendon driven continuum robot for use in minimal invasive surgery. Uses a semi-discretisation scheme to implicitly discretise in time and then solves the resulting ODEs.

Example input into dynamicShooting would be:

dt = 0.05; %Sets timestep size
t = 0:dt:5; %Sets time vector
Ns = 50; %Number of spacial segments
doPlot = "plot"; %Dynamic plotting

T_t = zeros(4,Nt); %Initialises tension vector
T_t(1,:) = 5*sin(5/30*pi()*t); %Sets one of the tensions to have a sinusoidal force input
Initial = "Straight"; %Input either "Straight" or "Static" for initial condition

dynamicShooting(T_t,Ns,doPlot, dt, t, Initial);
