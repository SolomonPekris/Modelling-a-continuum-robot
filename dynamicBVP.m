function out = dynamicShooting(T_t, N, Initial, dt, t, doPlot)
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
TRad = 8e-3; %Set Offset of Tendons from Backbone (m)
u_star = [0;0;0]; %Undeformed backbone linear speed configuration
v_star = [0;0;1]; %Undeformed backbone angular speed configuration
du_star = [0;0;0];
dv_star = [0;0;0];
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
xint = linspace(0,L,N);

% Set boundary conditions
p0 = [0;0;0];
q0 = [0;0;0];
w0 = [0;0;0];
h0 = [0.8256;0;0.8256;0];
if Initial == "Straight"
    %Initialises straight configurations
    %X and Z are general parameters, NOT coordinates
    % X = [p;R;n;m;q;w] and Z = [v;u]
     X = [linspace(0,L,N); zeros(18,N)];

      solution(:,:,1) = X(1:3,:);
      X = [p0;h0;zeros(6,1);q0;w0];
      Z = [0;0;1;0;0;0];
    
elseif Initial == "Static"
    T = T_t(:,1);
    X = staticShooting(T,N);
    Z = X(13:18,:);
    X = X(1:12,:);
    X = [X ; zeros(12,N)];
end

X_prev = X;
Z_prev = Z;

%BDF2 Coefficients. Change based on method used
c0=1.5/dt;
c1=-2/dt;
c2=0.5/dt;

%solution(:,:,1) = X(1:3,:);

sol.x = xint;
sol.y = X;
sol.yinit = X(1:19,1);

sol = bvpinit(xint,X(1:19,1));
j = 1;
visualise(doPlot) %Creates plot at current timestep
count = 0;
%Loops through time steps setting previous solution for n and m as
%the next guess
for j = 2:t_steps
    
    %Set history terms
    Xh = c1*X+c2*X_prev;
    Zh = c1*Z+c2*Z_prev;
    X_prev = X;
    Z_prev = Z;
   
    T = T_t(:,j); %Inputs tension at current timestep
    sol = bvp4c(@dynamicODE, @dynamicBVP, sol);
    count = count + 1;
    out = deval(sol,xint);
    
    X = out(1:19,1);
    h = X(4:7);
    h1=h(1);
    h2=h(2);
    h3=h(3);
    h4=h(4);
    R = eye(3) + 2/(h'*h) * ...
        [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
        h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
        h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];
    v0 = Kse\R.'*X(8:10) + v_star;
    u0 = Kbt\R.'*X(11:13) + u_star;
    Z(1:6) = [v0;u0];
    
%     for k = 1:N
%         h = X(4:7,k);
%         h1=h(1);
%         h2=h(2);
%         h3=h(3);
%         h4=h(4);
%         R = eye(3) + 2/(h'*h) * ...
%             [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
%             h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
%             h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];
%         
%         v0 = Kse\R.'*X(8:10,k) + v_star;
%         u0 = Kbt\R.'*X(11:13,k) + u_star;
%         Z(1:6,k) = [v0;u0];
%     end
    
    solution(1:3,:,j) = out(1:3,:);
    visualise(doPlot)

end
%solution = solution(:,:,t_steps);
%%
    function res = dynamicBVP(ra, rb)
       
        dpi = zeros(3,4);
        ni = zeros(3,4);
        mi = zeros(3,4);
           
        %Extract boundary variables
        nb = rb(8:10);
        mb = rb(11:13);

        h = rb(4:7);
     
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
        

        
        res = [ra(1:3,1)
            ra(4,1) - 0.8256
            ra(5,1)
            ra(6,1) - 0.8256
            ra(7,1)
            nb_error
            mb_error
            ra(14:19,1)];
        
    end

%%
    function [dgds] = dynamicODE(s,r)
  
        %Collect values from state variables
        h = r(4:7);
        n = r(8:10);
        m = r(11:13);
        q = r(14:16);
        w = r(17:19);
        h1=h(1);
        h2=h(2);
        h3=h(3);
        h4=h(4);
        R = eye(3) + 2/(h'*h) * ...
            [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
            h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
            h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];
        
        v = Kse\R.'*n + v_star;
        u = Kbt\R.'*m + u_star;
        
        Z = [v;u];

        %Time Derivatives
        Xt = c0*r + Xh;
        Zt = c0*Z + Zh;
        
        %Unpack derivatives
        vt = Zt(1:3,1);
        ut = Zt(4:6,1);
        qt = Xt(14:16,1);
        wt = Xt(17:19,1);
        
        % Intialise
        A = 0; B = 0; G = 0;
        H = 0; a = 0; b = 0;
        dp = zeros(3,4);
        
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

        c = rho*J*wt - cross(u,Kbt*(u - u_star))...
            - cross(v,Kse*(v-v_star)) - R.'*le - b + cross(w,rho*J*w);
        d = rho*Area*qt - cross(u,Kse*(v - v_star))...
            - R.'*fe - a + rho*Area*cross(w,q);
        
        dvdu = [Kse+A, G; B, Kbt + H]\[d; c];
        
        dvds = dvdu(1:3);
        duds = dvdu(4:6);
        ft = R*(a + A*dvds + G*duds);
        lt = R*(b + B*dvds + H*duds);
        %Calculates state vector dpds and dRds
        dpds = R*v;
        
        %Rearangment of eqn (27) from paper
        dqds = vt - cross(u,q) + cross(w,v);
        dwds = ut - cross(u,w);
        dnds = rho*Area*R*(cross(w,q) + qt) -fe-ft;
        dmds = R*(cross(w,rho*J*w) + rho*J*wt) -lt-le - cross(dpds,n);
        
        dhds =   [ 0, -u(1), -u(2), -u(3);
            u(1), 0, u(3), -u(2);
            u(2), -u(3), 0, u(1);
            u(3), u(2), -u(1), 0 ] * h/2 ;
  
        
        dgds = [dpds; dhds; dnds; dmds; dqds; dwds];
        
    end
%%
    function visualise(doPlot)
        if doPlot == "plot"
            plot(solution(1,:,j),solution(3,:,j)); title ( 'Cantilever Rod');
            xlabel('z (m)'); ylabel('x (m)');
            axis([0 1.1*L -0.8*L 0.8*L]);
            grid on; daspect([1 1 1]);
            legend(string(j*dt))
            drawnow;
        end
    end
%%
    function X = skew(x)
        %Function to calculate the skew of a vector
        X=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
        
    end
end
