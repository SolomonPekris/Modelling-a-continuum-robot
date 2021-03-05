dt = 0.05;
t = 0:dt:5;
Nt = length(t);
Ns = 40;
doPlot = "plot";
t_steps = max(t)/dt;

T_t = zeros(4,Nt);
5*sin(5/30*pi()*t)

Initial = "Straight";

solution = dynamicBVP(T_t, N, Initial, dt, t, doPlot);

x = solution(1,:,t_steps);
y = solution(2,:, t_steps);
z = solution(3,:,t_steps);

figure(1)
plot(x,z)
hold on
plot(out(1,:),out(3,:))%,'o','MarkerIndices',1:5:Ns)
xlabel('z')
ylabel('x')
% axis([0 1.1*0.24 -0.8*0.24 0.8*0.24]);
% daspect([1 1 1]);
for i = 1:1
    plot(Expi(i,1:3:21), Expi(i,3:3:21),'o')
end
hold off

grid on