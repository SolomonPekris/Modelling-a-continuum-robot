function [X, out] = staticShooting(T,Ns)
% Function that solves the governing equations based upon cosserat theory
% of rods and strings using shooting method

%Set Constants
options = optimset('Display','off');
L = 0.24; %Backbone Length (m)
Area = 5.02655e-7; %Backbone Cross Sectional Area (m^2)
TotalMass = 9e-3; %Total Continuum Robot Arm Mass (kg)
E = 168e9; %Modulus of Elasticity of Backbone (MPa)
Gs = 0.29*E; %Shear Modulus of Backbone (MPa)
g = -9.81; %Acceleration due to Gravity (m/s^2)
TRad = 8e-3; %Set Offset of Tendons from Backbone (m)
u_star = [0;0;0]; %Undeformed backbone linear speed configuration
v_star = [0;0;1]; %Undeformed backbone angular speed configuration
du_star = [0;0;0];
dv_star = [0;0;0];
fe = [0;0;TotalMass/L*g]; %Set Distributed Load on Backbone (Weight)
le = [0;0;0]; %Set Distributed Moment on Backbone
R0 = [0 0 1;0 1 0;-1 0 0]; %Rotational orientation of fixed end
h0 = [0.8256;0;0.8256;0];
p0 = [0;0;0]; %Position of fixed end in global coordinates
Ixx = 2.0106e-14; %Second moment of area
Iyy = Ixx;
Izz = 4.0212e-14; %Polar moment of inertia
Kse = diag([Gs*Area,Gs*Area,E*Area]); %Stiffness matrix shear and extension
Kbt = diag([E*Ixx,E*Iyy,E*Izz]); %Stiffness matrix= bending and torsion
s = linspace(0,L,Ns); % Creates spatial mesh
ri = [0 8e-3 0; -8e-3 0 0; 0 -8e-3 0; 8e-3 0 0];
ri = ri.';

%% Initialisation

%Int guess for n and m
rinita(1:6,:)=0.01;


%% Solve using shooting method

global X
fsolve(@staticShooting,rinita, options);
out = X(1:3,:);

%%
    function [dgds] = staticODE(s,r)
        % Int
        
        A = 0; B = 0; G = 0; H = 0; a = 0; b = 0; R = zeros(3,3); dp = zeros(3,4);
        
        %Collect values from state variables
        h = r(4:7);
        n = r(8:10);
        m = r(11:13);
        h1=h(1);
        h2=h(2);
        h3=h(3);
        h4=h(4);
        R = eye(3) + 2/(h'*h) * ...
            [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
            h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
            h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];
        %Calculate v and u from n and m
        v = Kse\R.'*n + v_star;
        u = Kbt\R.'*m + u_star;
        
        %Loop to calculate the variables for tendon load and moment
        for i = 1:4
            
            dp(:,i) = skew(u)*ri(:,i) + v;
            Ai = -T(i).*skew(dp(:,i))^2./norm(dp(:,i))^3;
            Bi = skew(ri(:,i))*Ai;
            A = A + Ai;
            B = B + Bi;
            G = G - Ai*skew(ri(:,i));
            H = H - Bi*skew(ri(:,i));
            a = a + Ai*(skew(u)*dp(:,i));
            b = b + skew(ri(:,i))*Ai*(skew(u)*dp(:,i));
            
        end
        
        
        
        %variables used to relate tendon loads and moments to dv and du
        c = Kbt*du_star - skew(u)*Kbt*(u - u_star) - skew(v)*Kse*(v-v_star) - R.'*le - b;
        d = Kse*dv_star - skew(u)*Kse*(v - v_star) - R.'*fe - a;
        
        dvdu = [Kse+A, G; B, Kbt + H]\[d; c];
        
        dv = [dvdu(1); dvdu(2); dvdu(3)];
        du = [dvdu(4); dvdu(5); dvdu(6)];
        
        %Distributed tendon load and moment
        ft = R*(a + A*dv + G*du);
        lt = R*(b + B*dv + H*du);
        
        %Calculates state variables
        dpds = R*v;
        dnds = -fe-ft;
        dmds = -lt-le - cross(dpds,n);
        dhds = [ 0  , -u(1), -u(2), -u(3);
            u(1),   0  ,  u(3), -u(2);
            u(2), -u(3),   0  ,  u(1);
            u(3),  u(2), -u(1),   0  ] * h/2;
        
        dgds = [dpds;dhds;dnds;dmds];
        
    end

%% Introduce System Boundary Conditions
    function res = staticShooting(guess)
        
        %Int
        dpi = zeros(3,4); ni = zeros(3,4); mi = zeros(3,4);
        
        %Guesses for n and m
        n0 = [guess(1);guess(2);guess(3)];
        m0 = [guess(4);guess(5);guess(6)];
        
        y0 = [p0
            h0
            n0
            m0];
        
        %Solve governing equations
        [~,X] = ode45(@staticODE,s,y0);
        X = X.';
        
        %Extract boundary variables
        nb = X(8:10, Ns);
        mb = X(11:13, Ns);
        
        h = X(4:7, Ns);
        h1 = h(1);
        h2 = h(2);
        h3 = h(3);
        h4 = h(4);
        Rb = eye(3) + 2/(h'*h) * ...
            [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
            h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
            h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];
        
        vb = Kse\Rb.'*nb + v_star;
        ub = Kbt\Rb.'*mb + u_star;
        % ri = tendonlocation(0.24,t_span);
        %Loop for boundary condition of tendon termination
        for i = 1:4
            dpi(:,i) = Rb*(skew(ub)*ri(:,i)+vb);
            
            ni(:,i) = T(i).*dpi(:,i)./norm(dpi(:,i));
            
            mi(:,i) = T(i)*skew(Rb*ri(:,i))*dpi(:,i)./norm(dpi(:,i));
        end
        
        %Error between actual BC and calculated BC
        nb_error = -nb - [sum(ni(1,:));sum(ni(2,:));sum(ni(3,:))];
        mb_error = -mb - [sum(mi(1,:));sum(mi(2,:));sum(mi(3,:))];
        
        %Error vector
        res = [nb_error; mb_error];
        
    end

%%
    function X = skew(x)
        % Function to calculate the skew of a vector
        X=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
        
    end
end
