close all
clear all
mu = 5;
G_target = 0.3; 
eta = 0.05;
kp = 0.3;
ki1 = 0.03;
ki2o = 0.45;
e_record = 0;
max_time = 100; %How long the simulation will run for.
dt = 0.01; %The time step we will use for our simple integration.
num = round(max_time/dt); % The number of time steps we are going to integrate for.
time = linspace(0,max_time,num+1); %A vector of times that will be used in plotting.

rng(1)

Brownian = true;
predictive = true;
PI_control = ~predictive;

[general_fig, alpha_est, n_est, G_est, freq, D_est, G_meas_fig] = deal(0,1,1,1,0,0,0);
[avg_filter, polyfit_filter, no_filter, sharp_edge_avg, preval_filter] = deal(0,0,0,0,1);
warning('off','all')
fig_time = 10;
fig_t = false;
ou_theta = 0.5;
ou_sigma = 0.3;
avg_filter_size = 100;  

x = 0.5;    % previous value weight

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
pulse = 50;
U_record(1,:)=pulse;
n_record(1,:)=n0;
D_record(1,:)=D0;
alpha_record(1,:)=alpha0;


dn_dt=@(n,D)(n.*(2*exp(-D)-1)-(n.^2)*exp(-D));
dD_dt=@(U,alpha,D)(U-alpha*D);
dalpha_dt=@(alpha,n,D)(alpha*n*sqrt(mu*D*alpha)*(1-exp(-D)));
dn_dt1=@(n,G)(2*G-n-n*G);


n_estimate=zeros(num+1,1);
D_estimate=zeros(num+1,1);
alpha_estimate=zeros(num+1,1);
G_measured=zeros(num+1,1);
G_filtered = zeros(num+1,1);


alpha_actual=zeros(num+1,1);
sigma=0.09;



sample = 0.1; %start of sampling
stop = 5.2; %end of sampling
ts = 0.1; %sampling interval
loop = ts/dt;
time_s=zeros(num+1,1);
i=1; j=2;                       % j = 2 to take initial values as given
n_estimate(1)=1;
G_measured(1)=1;
G_filtered(1)=1;


error = 0;


while alpha<1000
   n=n_record(i,:);
   D=D_record(i,:);
   alpha=alpha_record(i,:);
   G_curr=G(i,:);  

   nnew = n+dn_dt(n,D)*dt; 
   nnew(nnew<10^-10)=0;
  
  
   %U = (-D - log(G_target/(n)))/dt + alpha*D;
   if round(time(i)*100)/100 == round(sample*100)/100

      n_prev = n_estimate(j-1,:);
      G_prev = G_filtered(j-1,:);
      
      G_curr_m_W = normrnd(G_curr,sigma*G_curr); % noise
      
      if Brownian
        error = error + G_curr - G_curr_m_W; % brownian noise
        ou_error = -ou_theta*G_curr*ts + ou_sigma*error  ; 
      end
      
      
      G_measured(j,:) = G_curr_m_W + ou_error;
      G_measured(G_measured<10^-10,:)=0;
      G_curr_m = G_curr_m_W + ou_error;
      G_curr_m(G_curr_m<10^-10) = 0;
        
      % collecting previous samples
      
      k=0; y_fil = []; x_fil = [];
      if time(i) == 5.2            % debugging stop test time for breakpoints
          errorStop=true;
      end
      while (j-k) > 0 && k < avg_filter_size
         y_fil(k+1) = G_measured(j-k);                      %#ok<SAGROW>
         x_fil(k+1) = j-k;                                  %#ok<SAGROW>
         if  j-k-1>0 && sharp_edge_avg
            if ((G_measured(j-k-1)-G_measured(j-k))/G_measured(j-k)) > 0.1
                break
            end
         end
         k = k+1;
         
      end
      
      % moving average filter
      if avg_filter || sharp_edge_avg
          if j>2
              G_fil = mean(y_fil);
              G_filtered(j,:) = G_fil;
          else
              G_fil = G_curr_m;
              G_filtered(j,:) = G_fil;
          end
      elseif preval_filter
          G_fil = (1-x)*G_curr_m + x*G_prev;
          G_filtered(j,:) = G_fil; 
      elseif polyfit_filter    
          if j>2 
              G_fil = polyval(polyfit(x_fil, y_fil, 4), j);
              G_filtered(j,:) = G_fil;    
          else
              G_fil = G_curr_m;
              G_filtered(j,:) = G_fil;
          end
      elseif no_filter
          G_fil = G_curr_m;
          G_filtered(j,:) = G_fil;
      end
      % polyfit filter
      
      
      
      n_estimate(j,:) = n_prev + dn_dt1(n_prev,G_prev)*ts;
      n_now = n_estimate(j,:);
      
      D_estimate(j,:) = log(n_now/G_fil);
      D_estimate(D_estimate<10^-10) = 0;
      
      
      alpha_actual(j,:) = alpha;

      U_prev = U_record(i-1,:); 
      
      D_prev = D_estimate(j-1,:); D_now = D_estimate(j,:);
      
         
      alpha_estimate(j-1,:)=(U_prev*ts + D_prev - D_now)/(D_prev*ts);
      alpha_prev = alpha_estimate(j-1,:);
      
          
          
      if round(time(i)*100)/100 >= stop
         alpha_curr =  real(alpha_prev + dalpha_dt(alpha_prev,n_prev,D_prev)*ts);
         
         if predictive
            U_record(i,:) = max((((-log(n_now/G_fil) - log(G_target/(n_now)))/ts) + alpha_curr*log(n_now/G_fil)),0);
         elseif PI_control
             PI_error = G_fil - G_target;
             if abs(PI_error)>=0.15
                 ki2 = 0;
             else
                 ki2 = ki2o;
             end

             e_record = e_record + PI_error*ts;
             integral1=ki1*e_record;
             integral2=ki2*e_record;

             U_record(i,:) = (exp(eta*(kp*PI_error+integral1+integral2))-1)/eta;
             U_record(U_record<0)=0;
         end
      
         
      
      
      end
      
      time_s(j,:)=time(i);
      j=j+1; 
      sample = sample + ts; 
      
      
      
      
      
   elseif time(i)>dt*2
       U_record(i,:) = U_record(i-1,:);
   end
   
   U = real(U_record(i,:));
   alphanew = alpha+dalpha_dt(alpha,n,D)*dt;
   Dnew = D+dD_dt(U,alpha,D)*dt;
   Dnew(Dnew<10^-10)=0;
   
   %% Recording Final values each loop iteration
   n_record(i+1,:)=real(nnew);
   D_record(i+1,:)=real(Dnew);
   alpha_record(i+1,:)=real(alphanew);
   G(i+1,:)=nnew*exp(-Dnew);
   G(G<10^-10,:)=0;
i = i+1;   
end


if general_fig && ~fig_t
    figure
    subplot(2,2,1)
    plot(time(1:i-1),n_record(1:i-1));
    xlabel('Time');
    ylabel('n','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
    subplot(2,2,2)
    plot(time(1:i-1),real(U_record(1:i-1)));
    xlabel('Time');
    ylabel('U','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
    subplot(2,2,3)
    plot(time(1:i-1),D_record(1:i-1));
    xlabel('Time');
    ylabel('D','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
    subplot(2,2,4)
    plot(time(1:i-1),alpha_record(1:i-1));
    set(gca, 'YScale', 'log')
    xlabel('Time');
    ylabel('\alpha','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
elseif general_fig && fig_t
    figure
    subplot(2,2,1)
    plot(time(time<fig_time),n_record(time<fig_time));
    xlabel('Time');
    ylabel('n','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
    subplot(2,2,2)
    plot(time(time<fig_time),real(U_record(time<fig_time)));
    xlabel('Time');
    ylabel('U','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
    subplot(2,2,3)
    plot(time(time<fig_time),D_record(time<fig_time));
    xlabel('Time');
    ylabel('D','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
    subplot(2,2,4)
    plot(time(time<fig_time),alpha_record(time<fig_time));
    set(gca, 'YScale', 'log')
    xlabel('Time');
    ylabel('\alpha','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
end






if G_est 
    figure
    hold on
    if ~fig_t
        plot(time(1:i-1),G(1:i-1))
    else
        plot(time(time<fig_time),G(time<fig_time))
    end
    xlabel('Time');
    ylabel('G','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
    if ~fig_t
        plot(time_s(1:j-1),G_measured(1:j-1))
    else
        plot(time_s(time_s<fig_time),G_measured(time_s<fig_time))
    end
    if not(G_meas_fig)
        if ~fig_t
            plot(time_s(1:j-1),G_filtered(1:j-1))
        else
            plot(time_s(time_s<fig_time),G_filtered(time_s<fig_time))
        end
        leg = legend('G actual','G measured','G filtered');
    else
        leg = legend('G actual','G measured');
    end
    leg.ItemHitFcn = @hitcallback_ex1;
    ylim([0 1.5])
    hold off
end


if n_est
    figure
    hold on
    if ~fig_t
        plot(time(1:i-1),n_record(1:i-1));
    else
        plot(time(time<fig_time),n_record(time<fig_time));
    end
    xlabel('Time');
    ylabel('n','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
    if ~fig_t
        plot(time_s(1:j-1),n_estimate(1:j-1));
    else
        plot(time_s(time_s<fig_time),n_estimate(time_s<fig_time));
    end
    leg = legend('n actual','n estimate');
    leg.ItemHitFcn = @hitcallback_ex1;
    hold off
end

if D_est
    figure
    hold on
    if ~fig_t
        plot(time(1:i-1),D_record(1:i-1));
    else
        plot(time(time<fig_time),D_record(time<fig_time));
    end
    xlabel('Time');
    ylabel('D','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
    if ~fig_t
        plot(time_s(1:j-1),D_estimate(1:j-1));
    else
        plot(time_s(time_s<fig_time),D_estimate(time_s<fig_time));
    end
    leg = legend('D actual','D estimate');
    leg.ItemHitFcn = @hitcallback_ex1;
    hold off
end


if alpha_est
    figure
    hold on 
    if ~fig_t
        plot(time_s(1:j-2),alpha_actual(1:j-2));
    else
        plot(time_s(time_s<fig_time),alpha_actual(time_s<fig_time));
    end
    set(gca, 'YScale', 'log')
    xlabel('Time');
    ylabel('\alpha','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
    if ~fig_t
        plot(time_s(1:j-2),alpha_estimate(1:j-2));
    else
        plot(time_s(time_s<fig_time),alpha_estimate(time_s<fig_time));
    end
    leg = legend('\alpha actual','\alpha estimate');
    leg.ItemHitFcn = @hitcallback_ex1;
    hold off
end

if freq
    G_data = iddata([],G(2:i-1),0.01);
    GM_data = iddata([],G_measured(1:j-1),0.1);
    GF_data = iddata([],G_filtered(1:j-3),0.1);
    figure
    subplot(1,3,1)
    plot(fft(G_data))
    subplot(1,3,2)
    plot(fft(GM_data))
    subplot(1,3,3)
    plot(fft(GF_data))
end
