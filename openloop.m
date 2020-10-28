mu = 5;
max_time = 100; %How long the simulation will run for.
dt = 0.01; %The time step we will use for our simple integration.
n = round(max_time/dt); % The number of time steps we are going to integrate for.
time = linspace(0,max_time,n+1); %A vector of times that will be used in plotting.

n_record=zeros(n+1,1);
D_record=zeros(n+1,1);
alpha_record=zeros(n+1,1);
U_record=zeros(n+1,1);
U_record(time>10) = 0.1;
G = zeros(n+1,1);
%U_record = (exp(0.02*time)-1)';


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


for i=1:n 
   n=n_record(i,:);
   D=D_record(i,:);
   alpha=alpha_record(i,:);
   U = U_record(i,:);
   
   nnew = n+dn_dt(n,D)*dt; 
   Dnew = D+dD_dt(U,alpha,D)*dt;
   alphanew = alpha+dalpha_dt(alpha,n,D)*dt;
       
       
   
   %% Recording Final values each loop iteration
   n_record(i+1,:)=nnew;
   D_record(i+1,:)=Dnew;
   alpha_record(i+1,:)=alphanew;
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
