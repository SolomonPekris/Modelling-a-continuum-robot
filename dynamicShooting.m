function solution = dynamicShooting(T_t,N,doPlot, del_t, t,Initial)
%%Function that solves the governing equations based upon cosserat theory
%%of rods and strings using shooting method

t_steps = length(t);

%Set Constants
L = 0.24; %Backbone Length (m)
Area = 5.02655e-7; %Backbone Cross Sectional Area (m^2)
TotalMass = 9e-3; %Total Continuum Robot Arm Mass (kg)
E = 168e9; %Modulus of Elasticity of Backbone (MPa)
Gs = 0.29*E; %Shear Modulus of Backbone (MPa)
g = -9.81; %Acceleration due to Gravity (m/s^2)
u_star = [0;0;0]; %Undeformed backbone linear speed configuration
v_star = [0;0;1]; %Undeformed backbone angular speed configuration
fe = [0;0;TotalMass/L*g]; %Set Distributed Load on Backbone (Weight)
le = [0;0;0]; %Set Distributed Moment on Backbone
Ixx = 2.0106e-14; %Second moment of area
Iyy = Ixx;
Izz = 4.0212e-14; %Polar moment of inertia
Kse = diag([Gs*Area,Gs*Area,E*Area]); %Stiffness matrix shear and extension
Kbt = diag([E*Ixx,E*Iyy,Gs*Izz]); %Stiffness matrix= bending and torsion
ri = [0 8e-3 0; -8e-3 0 0; 0 -8e-3 0; 8e-3 0 0];
ri = ri.';
rho = 7860; %Density
J = diag([Ixx,Iyy,Izz]);
Bse = zeros(3); %Strain and extension damping
Bbt = zeros(3); %Bending and torsional damping
ds = L/(N-1);

% Set boundary conditions
p0 = [0;0;0];
h0 = [0.8256;0;0.8256;0];
q0 = [0;0;0];
w0 = [0;0;0];

guess(1:6,:) = 0.01; %Sets initial guess for n0 and m0

%Initialises first step as straight, or from given static position
if Initial == "Straight"
    %Initialize to straight configuration
    % y and z are general variables as in Eq (12)
    % y:=[p;h;n;m;q;w] and z:=[v;u]
    y = [linspace(0,L,N); zeros(18,N)];
    z = [zeros(2,N); ones(1,N); zeros(3,N)];
    
elseif Initial == "Static"
    T = T_t(:,1);
    myfun = @() staticShooting(T,N);
    y = staticShooting(T,N);
    guess = y(8:13,1);
    vinit = zeros(3,N);
    uinit = zeros(3,N);
    
    %Calculates v and u across for entire solution
    for i = 1:N
        h = y(4:7,i);
        h1=h(1);
        h2=h(2);
        h3=h(3);
        h4=h(4);
        R = eye(3) + 2/(h'*h) * ...
            [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
            h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
            h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];
        vinit(:,i) = Kse\R.'*y(8:10,i) + v_star;
        uinit(:,i) = Kbt\R.'*y(11:13,i) + u_star;
    end
    %Sets initial conditions
    z = [vinit;uinit];
    y = y(1:13,:);
    y = [y ; zeros(6,N)];
else
    error("Invalid initial condition")
end

y_prev = y;
z_prev = z;

%BDF2 Coefficients. Change based on method used
c0=1.5/del_t;
c1 =-2/del_t;
c2=0.5/del_t;
Kse_plus_c0_Bse_inv = (Kse+c0*Bse)^-1;
Kbt_plus_c0_Bbt_inv = (Kbt+c0*Bbt)^-1;
Kse_vstar = Kse*v_star;

solution(:,:,1) = y(1:3,:);
k = 0;
visualise(doPlot,k,del_t) %Creates plot at current timestep

%Loops through time steps setting previous solution for n and m as
%the next guess
for k = 2 : t_steps
    T = T_t(:,1);
    
    %Set history terms - Eq (5)
    yh = c1*y+c2*y_prev;
    zh = c1*z+c2*z_prev;
    y_prev = y;
    z_prev = z;
    
    %Midpoints are linearly interpolated for RK4
    yh_int = 0.5*(yh(:,1:end-1) + yh(:,2:end));
    zh_int = 0.5*(zh(:,1:end-1) + zh(:,2:end));
    
    %Shooting method solver call
    guess = fsolve(@getResidual, guess,options);
    
    solution(1:3,:,k) = y(1:3,:);
    
    visualise(doPlot,k,del_t);
end
%solution = solution(:,:,t_steps);
%%
    function E = getResidual(guess)
        %Reaction force and moment are guessed
        n0 = guess(1:3);
        m0 = guess(4:6);
        y(:,1) = [p0; h0; n0; m0; q0; w0];
        
        %Fourth-Order Runge-Kutta Integration
        for j = 1 : N-1
            yj = y(:,j);
            yhj_int = yh_int(:,j);
            [k1,z(:,j)]=ODE(yj,yh(:,j),zh(:,j));
            [k2,~]=ODE(yj+k1*ds/2,yhj_int,zh_int(:,j));
            [k3,~]=ODE(yj+k2*ds/2,yhj_int,zh_int(:,j));
            [k4,~]=ODE(yj+k3*ds,yh(:,j),zh(:,j+1));
            y(:,j+1) = yj + ds*(k1 + 2*(k2+k3) + k4)/6;
            %y(:,j+1) = yj + ds*k1; %Euler's Method
        end
        
        %Cantilever boundary conditions
        nb = y(8:10,N);
        mb = y(11:13,N);
        
        %Calculates rotation vector R from quaternion at tip
        h = y(4:7,N);
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
        
        %Loop for boundary condition of tendon termination
        for i = 1:4
            dpi(:,i) = Rb*(cross(ub,ri(:,i))+vb);
            
            ni(:,i) = T(i).*dpi(:,i)./norm(dpi(:,i));
            
            mi(:,i) = T(i)*skew(Rb*ri(:,i))*dpi(:,i)./norm(dpi(:,i));
        end
        
        %Error between actual BC and calculated BC
        nb_error = -nb - [sum(ni(1,:));sum(ni(2,:));sum(ni(3,:))];
        mb_error = -mb - [sum(mi(1,:));sum(mi(2,:));sum(mi(3,:))];
        
        %Sets residual
        E = [nb_error; mb_error];
    end

%%

    function [ys,z] = ODE(y,yh,zh)
        %Unpack y and z
        h = y(4:7);
        n = y(8:10);
        m = y(11:13);
        q = y(14:16);
        w = y(17:19);
        vh = zh(1:3);
        uh = zh(4:6);
        
        %Quaternion to Rotation
        h1=h(1);
        h2=h(2);
        h3=h(3);
        h4=h(4);
        R = eye(3) + 2/(h'*h) * ...
            [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
            h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
            h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];
        
        v=Kse_plus_c0_Bse_inv*(R'*n+Kse_vstar-Bse*vh);
        u=Kbt_plus_c0_Bbt_inv*(R'*m-Bbt*uh);
        z=[v;u];
        
        % Intialise
        A = 0; B = 0; G = 0;
        H = 0; a = 0; b = 0;
        dp = zeros(3,4);
        
        %Calculate and set time Derivatives
        yt = c0*y + yh;
        zt = c0*z + zh;
        
        vt = zt(1:3);
        ut = zt(4:6);
        qt = yt(14:16);
        wt = yt(17:19);
        
        for i = 1:4
            dp(:,i) = cross(u,ri(:,i)) + v;
            Ai = -T(i).*skew(dp(:,i))*skew(dp(:,i))./norm(dp(:,i))^3;
            Bi = skew(ri(:,i))*Ai;
            A = A + Ai;
            B = B + Bi;
            G = G - Ai*skew(ri(:,i));
            H = H - Bi*skew(ri(:,i));
            a = a + Ai*cross(u,dp(:,i));
            b = b + cross(ri(:,i),Ai*cross(u,dp(:,i)));
            
        end
        
        %Calculate constants c and d
        c = rho*J*wt - cross(u,Kbt*(u - u_star))...
            - cross(v,Kse*(v-v_star)) - R.'*le - b + cross(w,rho*J*w);
        d = rho*Area*qt - cross(u,Kse*(v - v_star))...
            - R.'*fe - a + rho*Area*cross(w,q);
        
        dvdu = [Kse+A, G; B, Kbt + H]\[d; c];
        
        dvds = dvdu(1:3);
        duds = dvdu(4:6);
        
        %Calculates total forces and moments acting along backbone
        ft = R*(a + A*dvds + G*duds);
        lt = R*(b + B*dvds + H*duds);
        
        %Calculates state vector dpds and dRds
        dpds = R*v;
        
        %Solves ODEs        
        dqds = vt - cross(u,q) + cross(w,v);
        dwds = ut - cross(u,w);
        dnds = rho*Area*R*(cross(w,q) + qt) -fe-ft;
        dmds = R*(cross(w,rho*J*w) + rho*J*wt) -lt-le - cross(dpds,n);
        
        dhds = [ 0  , -u(1), -u(2), -u(3);
            u(1),   0  ,  u(3), -u(2);
            u(2), -u(3),   0  ,  u(1);
            u(3),  u(2), -u(1),   0  ] * h/2;
        
        ys = [dpds;dhds;dnds;dmds;dqds;dwds];
        
    end
%%
    function visualise(doPlot,k,del_t)
        %Plots dynamic robot position in 2D
        if doPlot == "plot"
            plot(y(1,:),y(3,:)); title ( 'Cantilever Rod');
            xlabel('z (m)'); ylabel('x (m)');
            axis([0 1.1*L -0.8*L 0.8*L]);
            grid on; daspect([1 1 1]);
            legend(string(k*del_t))
            drawnow;
        end
    end
%%
    function X = skew(x)
        %Function to calculate the skew of a vector
        X=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
        
    end
end
