clc;
clear;
step=5000;
Volume_allot=2.4*10^8;
Load=[450 450 500 500 550 550 650 650 700 700 800 800 850 850 800 800 750 750 700 700 600 600 500 500  ];

Ph_max=700;
Q_max=(30*Ph_max+0.02*Ph_max*Ph_max);
Volume(25)=0;
Vmax(25)=0;
Vmin(25)=0;
for hour=24:-1:1 %%upper limit detection
    Vmax(hour)=Vmax(hour+1)+Q_max*60*60;
    if Vmax(hour)>Volume_allot
        Vmax(hour)=Volume_allot;
    end
end
Vmin(1)=Volume_allot;
for hour=2:1:24%%lower limit detect
    Vmin(hour)=Vmin(hour-1)-Q_max*60*60;
    if Vmin(hour)<0
        Vmin(hour)=0;
    end
end
Q_initial=(Volume_allot/24);
Volume(25)=0;
Volume_above(25)=0;
Volume_below(25)=0;
for hour=24:-1:1  %%avg/mid line detection
    Volume(hour)=Volume(hour+1)+Q_initial;
    Volume_below(hour)=Volume(hour)-step;
    Volume_above(hour)=Volume(hour)+step;
end
Volume_initial=Volume;
 hour=1:25;
 plot(hour,Vmax,'b',hour,Vmin,'r',hour,Volume,'g');
axis([1 25 1 600000000]);
 grid on;
%%
%%%combine thermal units into 1;
gen_data=[ 0.004 12.2 300 50 500
           0.0035 12.3 200 40 400
          ]; 
%%represent the given thermal unit by an equivalent unit
%Calculating values of lambda,there will be 2 lambda for each units%
     count=2*2; 
     lamda=zeros(count,1);
     power_fordifflamdas=zeros(count,1);
     j=1;
     a=1;
     while(j<=count)
         lamda(j)=2*gen_data(a,1)*gen_data(a,4)+gen_data(a,2);
         lamda(j+1)=2*gen_data(a,1)*gen_data(a,5)+gen_data(a,2);
         a=a+1;
         j=j+2;
     end
     a=1;
     while(a<=count)
     for j=1:1:2 
         powerofunit(j)=(lamda(a)-gen_data(j,2))/(2*gen_data(j,1));
         if  powerofunit(j)>gen_data(j,5)
              powerofunit(j)=gen_data(j,5);
         end
             if  powerofunit(j)<gen_data(j,4)
              powerofunit(j)=gen_data(j,4);
             end
             power_fordifflamdas(a)=power_fordifflamdas(a)+powerofunit(j);
     end
         a=a+1;
     end
     lamda=sort(lamda);
     power_fordifflamdas=sort(power_fordifflamdas);
     t=table(power_fordifflamdas,lamda);
     figure(1);
     plot(power_fordifflamdas,lamda);
     xlabel('Total Power');
     ylabel('System Lamda');
     grid on;
     %%
     %%%find initial G
    % G(1:24)=20;
for hour=1:1:24
    A=0.02; %%from given Q equation
    B=30;
    C=-Q_initial;
init_ph= (-B-sqrt(abs(B*B-4*A*C)))/(2*A); 
        a=0.0002;
            c=0.0004*init_ph*init_ph+Load(hour)-1.0006*init_ph;
            %used loss equation and basic power consumtion equation to find
            %the 2nd order eqn and using its root deriving formula to find
            %pt
            b=0.0003*init_ph-1.0005;%chg%chg
            ini_Ps(hour)=(-b-sqrt(abs(b*b-4*a*c)))/(2*a);
             PF_hydro__ini(hour)=1/(1-2*0.0004*init_ph+0.0003*ini_Ps(hour)-0.0006);
            PF_steam_ini(hour)=1/(1-2*0.0002*ini_Ps(hour)+0.0003*init_ph-0.0005);
           G(hour)=interp1(power_fordifflamdas,lamda,ini_Ps(hour))*(PF_steam_ini(hour)/PF_hydro__ini(hour));       
end
Power_hydro_old(1:24)=init_ph;
  %%
Q=zeros(3,3,24);
Ph=zeros(3,3,24);
Ps=zeros(3,3,24);
cost=zeros(3,3,24);
Ploss=zeros(3,3,24);
iteration=2000;
for loop_no=1:1:iteration
  for hour=24:-1:1
  if hour==24
      for state_ini=1:1:3
            if state_ini==1
                Volume_first=Volume_above;
            end
            if state_ini==2
                Volume_first=Volume;
            end
            if state_ini==3
                Volume_first=Volume_below;
            end
           state_fin=2;
            if state_fin==1
                Volume_second=Volume_above;
            end
            if state_fin==2
                Volume_second=Volume;
            end
            if state_fin==3
                Volume_second=Volume_below;
            end
            Q(state_ini,state_fin,hour)=(Volume_first(hour)-Volume_second(hour+1))/3600;
             %%%from Q equation
             
  
            C=-Q(state_ini,state_fin,hour); 
            A=0.02;B=30;
            Ph1=(-B+sqrt(abs(B*B-4*A*C)))/(2*A);
            Ph(state_ini,state_fin,hour)=Ph1; %got ph now need ps
             %used loss equation and basic power consumtion equation to find
            %the 2nd order eqn and using its root deriving formula to find
            %pt
            a=0.0002;
            c=0.0004*Ph(state_ini,state_fin,hour)*Ph(state_ini,state_fin,hour)+Load(hour)-1.0006*Ph(state_ini,state_fin,hour);%chang
            b=0.0003*Ph(state_ini,state_fin,hour)-1.0005;%chang
            Ps1=(-b-sqrt(abs(b*b-4*a*c)))/(2*a);
            Ps(state_ini,state_fin,hour)=Ps1;
            % ploss from quesstion 
            Ploss(state_ini,state_fin,hour)=0.0004*Ph(state_ini,state_fin,hour)*Ph(state_ini,state_fin,hour)+0.0002*Ps(state_ini,state_fin,hour)*Ps(state_ini,state_fin,hour)+0.0003*Ph(state_ini,state_fin,hour)*Ps(state_ini,state_fin,hour)-0.0006*Ph(state_ini,state_fin,hour);
            Max_cost(state_ini,hour)=G(hour)*Ph(state_ini,state_fin,hour);%chang  3 stage so 3 cost each hr
      end     
  end
  if hour>=2&&hour<=23
      for state_ini=1:1:3
            if state_ini==1
                Volume_first=Volume_above;
            end
            if state_ini==2
                Volume_first=Volume;
            end
            if state_ini==3
                Volume_first=Volume_below;
            end
           for state_fin=1:1:3
            if state_fin==1
                Volume_second=Volume_above;
            end
            if state_fin==2
                Volume_second=Volume;
            end
            if state_fin==3
                Volume_second=Volume_below;
            end
            Q(state_ini,state_fin,hour)=(Volume_first(hour)-Volume_second(hour+1))/3600;
              
              
                 C=-Q(state_ini,state_fin,hour);
            A=0.02;B=30;
            
            Ph1=(-B+sqrt(abs(B*B-4*A*C)))/(2*A);
            Ph(state_ini,state_fin,hour)=Ph1;
              
           %used loss equation and basic power consumtion equation to find
            %the 2nd order eqn and using its root deriving formula to find
            %pt
            a=0.0002;%chang
            c=0.0004*Ph(state_ini,state_fin,hour)*Ph(state_ini,state_fin,hour)+Load(hour)-1.0006*Ph(state_ini,state_fin,hour);%chang
            b=0.0003*Ph(state_ini,state_fin,hour)-1.0005;%
     
            Ps1=(-b-sqrt(abs(b*b-4*a*c)))/(2*a);
            Ps(state_ini,state_fin,hour)=Ps1;
            Ploss(state_ini,state_fin,hour)=0.0004*Ph(state_ini,state_fin,hour)*Ph(state_ini,state_fin,hour)+0.0002*Ps(state_ini,state_fin,hour)*Ps(state_ini,state_fin,hour)+0.0003*Ph(state_ini,state_fin,hour)*Ps(state_ini,state_fin,hour)-0.0006*Ph(state_ini,state_fin,hour);

            cost(state_ini,state_fin,hour)=G(hour)*Ph(state_ini,state_fin,hour)+Max_cost(state_fin,hour+1);
           end
           [Max_cost(state_ini,hour) Next_state(state_ini,hour)]=max(cost(state_ini,:,hour));
      end
  end
  if hour==1
      for state_fin=1:1:3
           if state_fin==1
                Volume_second=Volume_above;
            end
            if state_fin==2
                Volume_second=Volume;
            end
            if state_fin==3
                Volume_second=Volume_below;
            end
            Volume_first=Volume;
            Qist(hour,state_fin)=(Volume_first(hour)-Volume_second(hour+1))/3600;
            
             
            Cist=-Q(hour,state_fin); %chng
            Aist=0.02;Bist=30;%chng A=0.02;B=30;
            Ph1ist=(-Bist+sqrt(abs(Bist*Bist-4*Aist*Cist)))/(2*Aist);
            Phist(hour,state_fin)=Ph1ist;
            %used loss equation and basic power consumtion equation to find
            %the 2nd order eqn and using its root deriving formula to find
            %pt
            aist=0.0002;
            cist=0.0004*Phist(hour,state_fin)*Phist(hour,state_fin)+Load(hour)-1.0006*Phist(hour,state_fin);
            bist=0.0003*Phist(hour,state_fin)-1.0005;
            Ps1ist=(-bist-sqrt(abs(bist*bist-4*aist*cist)))/(2*aist);
            Psist(hour,state_fin)=Ps1ist;
            Plossist(hour,state_fin)=0.0004*Phist(hour,state_fin)*Phist(hour,state_fin)+0.0002*Psist(hour,state_fin)*Psist(hour,state_fin)+0.0003*Phist(hour,state_fin)*Psist(hour,state_fin)-0.0006*Phist(hour,state_fin);
            costist(hour,state_fin)=G(hour)*Phist(hour,state_fin)+Max_cost(state_fin,hour+1);
      end
  end
  end
[Value position]=max(costist);
new_state(1)=2;
new_state(2)=position;
new_state(25)=2;
for i=3:1:24
    new_state(i)=Next_state(new_state(i-1),i-1);
end
Power_hydro(1)=Phist(new_state(2));
Power_steam(1)=Psist(new_state(2));
Power_loss(1)=Plossist(new_state(2));
for i=2:1:23
    Power_hydro(i)=Ph(new_state(i),new_state(i+1),i);
    Power_steam(i)=Ps(new_state(i),new_state(i+1),i);
    Power_loss(i)=Ploss(new_state(i),new_state(i+1),i);
end
Power_hydro(24)=Ph(new_state(i),2,24);
Power_steam(24)=Ps(new_state(i),2,24);
Power_loss(24)=Ploss(new_state(i),2,24);
for hour=1:1:24
    
    PF_hydro(hour)=1/(1-2*0.0004*Power_hydro(hour)-0.0003*Power_steam(hour)-0.0006); % usinf d/dph *ploss
    PF_steam(hour)=1/(1-2*0.0002*Power_steam(hour)-0.0003*Power_hydro(hour)-0.0005);% usinf d/dps *ploss
    G(hour)=interp1(power_fordifflamdas,lamda,Power_steam(hour))*( PF_steam(hour)/PF_hydro(hour));
end
Volume(1)=Volume(1);
for hour=2:1:24
    if new_state(hour)==1
        Volume(hour)=Volume(hour)+step;
    end
    if new_state(hour)==2
        Volume(hour)=Volume(hour);
    end
    if new_state(hour)==3
        Volume(hour)=Volume(hour)-step;
    end
    Volume_above(hour)=Volume(hour)+step;
    Volume_below(hour)=Volume(hour)-step;
end
Volume(25)=0;
Volume_above(25)=0;
Volume_below(25)=0;
error=max(abs((Power_hydro-Power_hydro_old)));
end
hour=1:25;
figure(2);
plot(hour,Vmax,hour,Vmin,hour,Volume,hour,Volume_initial);
axis([1 25 1 700000000]);
xlabel('Hour');
ylabel('Volume');
legend('Upper Limit','Lower Limit','Voptimum','Vinitial');
grid on;
%%
%%%find indivudal power outputs
for hour=1:1:24
lamda_sys(hour)=interp1(power_fordifflamdas,lamda,Power_steam(hour));
P1(hour)=(lamda_sys(hour)-gen_data(1,2))/(2*gen_data(1,1));
P2(hour)=(lamda_sys(hour)-gen_data(2,2))/(2*gen_data(2,1));
 PF_hydro(hour)=1/(1-2*0.0004*Power_hydro(hour)-0.0003*Power_steam(hour)-0.0006);% usinf d/dph *ploss
    PF_steam(hour)=1/(1-2*0.0002*Power_steam(hour)-0.0003*Power_hydro(hour)-0.0005);% usinf d/dps *ploss
    

cost1(hour)=gen_data(1,1)*P1(hour)*P1(hour)+gen_data(1,2)*P1(hour)+gen_data(1,3);
cost2(hour)=gen_data(2,1)*P2(hour)*P2(hour)+gen_data(2,2)*P2(hour)+gen_data(2,3);

Total_cost(hour)=cost1(hour)+cost2(hour);

    Discharge(hour)=0.02*(Power_hydro(hour))*(Power_hydro(hour))+30*(Power_hydro(hour))*60*60;

if hour==1
    Volume_state(hour)=Volume_allot;
else   
Volume_state(hour)=Volume_state(hour-1)-Discharge(hour-1);
end
end


