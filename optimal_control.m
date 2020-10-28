close all
mu = 5;
G_target = 0.7; 
eta = 0.05;
kp = 0.3;
ki1 = 0.03;
ki2o = 0.45;
max_time = 100; %How long the simulation will run for.
dt = 0.01; %The time step we will use for our simple integration.
num = round(max_time/dt); % The number of time steps we are going to integrate for.
time = linspace(0,max_time,num+1); %A vector of times that will be used in plotting.

n_record=zeros(num+1,1);
D_record=zeros(num+1,1);
alpha_record=zeros(num+1,1);
U_record=zeros(num+1,1);
G = zeros(num+1,1);
e_record = 0;
n0=1;
D0=0;
alpha0=0.1;

n_record(1,:)=n0;
D_record(1,:)=D0;
alpha_record(1,:)=alpha0;

U_test = 0:0.01:100;

dn_dt=@(n,D)(n.*(2*exp(-D)-1)-(n.^2)*exp(-D));
dD_dt=@(U,alpha,D)(U-alpha*D);
dalpha_dt=@(alpha,n,D)(alpha*n*sqrt(mu*D*alpha)*(1-exp(-D)));

G(1,:)=n0*exp(-D0);


for i=1:num 
   n=n_record(i,:);
   D=D_record(i,:);
   alpha=alpha_record(i,:);
   G_curr=G(i,:);

   nnew = n+dn_dt(n,D)*dt; 
   nnew(nnew<0)=0;
   D_test = D+dD_dt(U_test,alpha,D)*dt;
   D_test(D_test<0)=0;
   alphanew = alpha+dalpha_dt(alpha,n,D)*dt;

   G_test = nnew*exp(-D_test);
   [~, ind] = min(abs(G_test-G_target));
   U = U_test(ind);
   Dnew = D+dD_dt(U,alpha,D)*dt;
   Dnew(Dnew<0)=0;
   
   %% Recording Final values each loop iteration
   n_record(i+1,:)=nnew;
   D_record(i+1,:)=Dnew;
   alpha_record(i+1,:)=alphanew;
   U_record(i+1,:)=U;
   G(i+1,:)=nnew*exp(-Dnew);
end

figure
subplot(2,2,1)
plot(time,n_record);
xlabel('Time');
ylabel('n','Rotation',0);
subplot(2,2,2)
plot(time,U_record);
xlabel('Time');
ylabel('U','Rotation',0);
subplot(2,2,3)
plot(time,D_record);
xlabel('Time');
ylabel('D','Rotation',0);
subplot(2,2,4)
plot(time,alpha_record);
set(gca, 'YScale', 'log')
xlabel('Time');
ylabel('\alpha','Rotation',0);
figure
plot(time,G)
xlabel('Time');
ylabel('G','Rotation',0);