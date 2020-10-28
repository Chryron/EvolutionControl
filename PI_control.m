mu = 5;
G_target = 0.7; 
eta = 0.0035;
kp = 0.01;
ki1 = 0.01;
ki2o = 0.35;
max_time = 100; %How long the simulation will run for.
dt = 0.01; %The time step we will use for our simple integration.
num = round(max_time/dt); % The number of time steps we are going to integrate for.
time = linspace(0,max_time,num+1); %A vector of times that will be used in plotting.

n_record=zeros(num+1,1);
D_record=zeros(num+1,1);
alpha_record=zeros(num+1,1);
U_record=zeros(num+1,1);
G = zeros(num+1,1);


e_record = zeros(num+1,1);

n0=1;
D0=0;
alpha0=0.1;

n_record(1,:)=n0;
D_record(1,:)=D0;
alpha_record(1,:)=alpha0;


dn_dt=@(n,D)(n.*(2*exp(-D)-1)-(n.^2)*exp(-D));
dD_dt=@(U,alpha,D)(U-alpha*D);
dalpha_dt=@(alpha,n,D)(alpha*n*sqrt(mu*D*alpha)*(1-exp(-D)));

G(1,:)=n0*exp(-D0);


for i=1:num 
   n=n_record(i,:);
   D=D_record(i,:);
   alpha=alpha_record(i,:);
   G_curr=G(i,:);
   Gmax = max(G(1:i));
   if eta>=0.15*Gmax
       ki2 = 0;
   else
       ki2 = ki2o;
   end
   
   e_record(i) = G_target - G_curr; 
   integral1=integral1+dt*ki1*error;
   integral1=integral2+dt*ki2*error*(abs(error-target)<0.2);
   
   U = (exp(eta*(kp*e_record(i)+ki1*trapz(e_record(1:i))+ki2*trapz(e_record(1:i))))-1)/eta;
   U(U<0)=0;
   nnew = n+dn_dt(n,D)*dt; 
   nnew(nnew<0)=0;
   Dnew = D+dD_dt(U,alpha,D)*dt;
   Dnew(Dnew<0)=0;
   alphanew = alpha+dalpha_dt(alpha,n,D)*dt;
   
       
   
   %% Recording Final values each loop iteration
   n_record(i+1,:)=nnew;
   D_record(i+1,:)=Dnew;
   alpha_record(i+1,:)=alphanew;
   U_record(i+1,:)=U;
   G(i+1,:)=dn_dt(nnew,Dnew);
end

figure
subplot(2,2,1)
plot(time,n_record);
xlabel('Time');
ylabel('n');
subplot(2,2,2)
plot(time,U_record);
xlabel('Time');
ylabel('U');
subplot(2,2,3)
plot(time,D_record);
xlabel('Time');
ylabel('D');
subplot(2,2,4)
plot(time,alpha_record);
set(gca, 'YScale', 'log')
xlabel('Time');
ylabel('alpha');
figure
plot(time,G)
xlabel('Time');
ylabel('G');