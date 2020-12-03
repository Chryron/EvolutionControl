close all
mu = 5;
G_target = 0.3; 
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

n0=1;
D0=0;
alpha0=0.1;
alpha = alpha0;

G(1,:)=n0;
U_record(1,:)=50;
n_record(1,:)=n0;
D_record(1,:)=D0;
alpha_record(1,:)=alpha0;


dn_dt=@(n,D)(n.*(2*exp(-D)-1)-(n.^2)*exp(-D));
dD_dt=@(U,alpha,D)(U-alpha*D);
dalpha_dt=@(alpha,n,D)(alpha*n*sqrt(mu*D*alpha)*(1-exp(-D)));
dn_dt1=@(n,G)(2*G-n-n*G);


n_measured=zeros(num+1,1);
D_estimate=zeros(num+1,1);
alpha_estimate=zeros(num+1,1);
G_measured=zeros(num+1,1);



alpha_actual=zeros(num+1,1);
sigma=0.01;



sample = 0.1; %start of sampling
stop = 5.1; %end of sampling
ts = 0.1; %sampling interval
loop = ts/dt;
time_s=zeros(num+1,1);
i=1; j=1;

while alpha<10^2
   n=n_record(i,:);
   D=D_record(i,:);
   alpha=alpha_record(i,:);
   G_curr=G(i,:);  

   nnew = n+dn_dt(n,D)*dt; 
   nnew(nnew<0)=0;
  
  
   %U = (-D - log(G_target/(n)))/dt + alpha*D;
   if round(time(i)*100)/100 == round(sample*100)/100
      if j==1
          n_prev = 1;
          G_prev = 1;
      else
          n_prev = n_measured(j-1,:);
          G_prev = G_measured(j-1,:);
      end
      G_measured(j,:) = normrnd(G_curr,sigma); % noise
      n_measured(j,:) = n_prev + dn_dt1(n_prev,G_prev)*ts;
      n_now = n_measured(j,:);
      D_estimate(j,:) = log(n_now/G_curr);
      %alpha_actual(j,:) = alpha;
      if j>1
          U_prev = U_record(i-1,:); 
          D_prev = D_estimate(j-1,:); D_now = D_estimate(j,:);
          alpha_estimate(j-1,:)=(U_prev*ts + D_prev - D_now)/(D_prev*ts);
          alpha_prev = alpha_estimate(j-1,:);
      end
      
      if round(time(i)*100)/100 >= stop
         alpha_curr =  alpha_prev + dalpha_dt(alpha_prev,n_prev,D_prev)*ts;
         U_record(i,:) = (-D_now - log(G_target/(n_now)))/ts + alpha_curr*D_now;
      end
      time_s(j,:)=time(i);
      j=j+1; 
      sample = sample + ts; 
   elseif time(i)>dt*2
       U_record(i,:) = U_record(i-1,:);
   end
   
   U = U_record(i,:);
   alphanew = alpha+dalpha_dt(alpha,n,D)*dt;
   Dnew = D+dD_dt(U,alpha,D)*dt;
   Dnew(Dnew<0)=0;
   
   %% Recording Final values each loop iteration
   n_record(i+1,:)=nnew;
   D_record(i+1,:)=Dnew;
   alpha_record(i+1,:)=alphanew;
   G(i+1,:)=nnew*exp(-Dnew);

i = i+1;   
end

figure
subplot(2,2,1)
plot(time(1:i-1),n_record(1:i-1));
xlabel('Time');
ylabel('n','Rotation',0);
subplot(2,2,2)
plot(time(1:i-1),U_record(1:i-1));
xlabel('Time');
ylabel('U','Rotation',0);
subplot(2,2,3)
plot(time(1:i-1),D_record(1:i-1));
xlabel('Time');
ylabel('D','Rotation',0);
subplot(2,2,4)
plot(time(1:i-1),alpha_record(1:i-1));
set(gca, 'YScale', 'log')
xlabel('Time');
ylabel('\alpha','Rotation',0);

figure
subplot(1,2,1)
plot(time(1:i-1),G(1:i-1))
xlabel('Time');
ylabel('G','Rotation',0);
subplot(1,2,2)
plot(time_s(1:j-1),G_measured(1:j-1))
ylim([0.15,1])
xlabel('Time');
ylabel('G measured','Rotation',0);



figure
subplot(1,2,1)
plot(time_s(1:j-1),n_measured(1:j-1));
xlabel('Time');
ylabel('n est','Rotation',0);
subplot(1,2,2)
plot(time(1:i-1),n_record(1:i-1));
xlabel('Time');
ylabel('n act','Rotation',0);
% figure
% subplot(2,2,1)
% plot(time_s(1:j-1),n_measured(1:j-1));
% xlabel('Time');
% ylabel('n est','Rotation',0);
% subplot(2,2,2)
% plot(time_s(1:j-1),G_measured(1:j-1));
% xlabel('Time');
% ylabel('G est','Rotation',0);
% subplot(2,2,3)
% plot(time_s(1:j-1),D_estimate(1:j-1));
% xlabel('Time');
% ylabel('D est','Rotation',0);
% subplot(2,2,4)
% plot(time_s(1:j-2),alpha_estimate(1:j-2));
% set(gca, 'YScale', 'log')
% xlabel('Time');
% ylabel('\alpha est','Rotation',0);
% 
% figure
% subplot(1,2,1)
% plot(time_s(1:j-2),alpha_estimate(1:j-2));
% set(gca, 'YScale', 'log')
% xlabel('Time');
% ylabel('\alpha est','Rotation',0);
% subplot(1,2,2)
% plot(time_s(1:j-2),alpha_actual(1:j-2));
% set(gca, 'YScale', 'log')
% xlabel('Time');
% ylabel('\alpha act','Rotation',0);
