
%
% make a read file for LAMMPS

clear; % clear all variables, RAM etc
clf; % clear (current) figure
tic
%numbers

Nmon1=2 %36;
Nmon2=64-Nmon1 %35;
Nmon3=0;


q=1;

Lz=290;

Nbrush = 1; %% TWO BRUSHES OR 1 ???

solventtype=0; %% no solvent (0), monomeric (1) or dimeric solvent (2)

if ((solventtype ~=1) && (solventtype ~= 0) && (solventtype ~= 2)) 
    solventtype
    fprintf('Choose value 0 (none), 1 (monomeric) or 2 (dimeric) for solventtype!\n'); 
    return;
end   

%%% gut: "factor" ist quadrat einer natuerlichen zahl oder
%%% eines bruchs natürlicher zahlen
factor=1; %4   %16/9; %25/9; %4; %% increase initial density by this factor

%Nmontot=Nmon1+Nmon2*2+Nmon3*4;
Nmontot=Nmon1+Nmon2*q%+Nmon3*4
%return

%return
%sigma=0.075 %0.0;5 %grafting density
sigma_final=0.100 %0.025 %2
sigma=sigma_final/factor

if ((Nbrush ~=1) && (Nbrush ~= 2)) 
    Nbrush
    fprintf('Choose value 1 or 2 for Nbrush!\n'); 
    return;
end   

if (Nmon1==0)
    fprintf('At least one monomer of stem polymer needed! (Nmon1 ~= 0 necessary)\n'); 
    return;
end    

if (((Nmon1 == 0) && (Nmon2 ~= 0)) | ((Nmon1 == 0) && (Nmon3 ~= 0)) | ((Nmon2 == 0) && (Nmon3 ~= 0)))
    fprintf('Lower generations of branching should exist before generating higher generations!\n'); 
    return;
end

rho=0.85;
Lz_final=12.0; %100.00  %12.0%17.5;

Lz_final1= 12.0
Lz_final2= 14.25
Lz_final3= 17.5

%Lx_final=63;
%Ly_final=43.6476;

%Lx_init=70;
%Ly_init=50;


%% number of rows with atoms in respective dimension
%xrows=50; %120
%yrows=10; %24

%xrows=180;
%yrows=40;
%xrows=152*3/4;
%yrows=32*3/4;

%xrows=102;
%yrows=36;

%xrows=118;
%yrows=27;

%xrows=152*3/4-10;
%yrows=32*3/4+5;
xrows=(152*3/4-14)*1.5
yrows=(32*3/4+6)*1.5

%return
 %xrows=30;
 %yrows=10;
nnn=xrows*yrows; %% number of lattice atoms

%% array to save coordinates of l21;attice 
b = struct('x',zeros(nnn,1),'y',zeros(nnn,1),'z',zeros(nnn,1));

%% increment of lattice 
dx=0.52500000;
dy=0.9093250000;
%dy=0.9093;

%% starting from zero
xstart=0;
ystart=0;
zstart=0;

%% calculating coordinates for lattice atoms
index=1;
for i=1:xrows
    for j=1:yrows
        b.x(index)=xstart+(i-1)*dx;
        b.y(index)=ystart+(j-1)*2*dy+mod(i+1,2)*dy;
        b.z(index)=zstart;
        index=index+1;
    end
end

%% make it symmetrical around zero in x and y
b.x(:)=b.x(:)-xrows*dx/2;
b.y(:)=b.y(:)-yrows*dy;

%% determine lattice size in x and y
xmin=min(b.x);
xmax=max(b.x);

ymin=min(b.y);
ymax=max(b.y);

%% calculate box size (adding increment so, that it can be used in period boundary conditions)
% xbox=xmax-xmin+dx
% ybox=ymax-ymin+dy
% nnn

Lx=xmax-xmin+dx
Ly=ymax-ymin+dy
Nc = floor(sigma*Lx*Ly*factor)+1

 Nc = Nc -1 + 0
%Nc = floor(sigma*Lx_final*Ly_final*factor)+1
%Nc=2
graftingpoints = zeros(2*Nc,2);
graftingcoords = zeros(2*Nc,2);



%% if we have Nc (number of chains in total)
%% we can sample from the Schulz-Zimm (SZ) distribution
%% to get a polydisperse distribution of chains

%% first, define the parameters of the SZ-d.

K = 50
polyd = 1 + 1/K;
Nn = Nmon2+Nmon1 % e.g. 128


%% plot the function
%my_x = [0:0.1:300]
my_x = [0:1:15*Nn];

%% Schulz-Zimm  - there is a problem with the gamma function here
SZ =   K.^K  * my_x.^(K-1) / (gamma(K) * Nn.^(K))  .* exp(-1.0 * K * my_x / Nn);
%SZ_0 = K.^K  * my_x.^(K-1) / (Nn.^(K))  .* exp(-1.0 * K * my_x / Nn)
%SZ_0 = K.^K  * my_x.^(K-1) / (gamma(K) * Nn.^(K))
%SZ_1 = exp(-1.0 * K * my_x / Nn)

%my_gammaln = gammaln(K)



sum(SZ);

%SZ_c = round(SZ*Nc*0) % discretized SZ distribution with integral = Nc

%Nc = 60000

lll = 0
while (lll <  10000)

%SZ_c = round(SZ*Nc*0.15); % discretized SZ distribution with integral = Nc 
SZ_c = zeros(1,length(my_x));




%% lets get this by pulling random numbers according to the SZ distr. above
%SZ_c_rd = zeros(1,2*Nn);

%length(my_x)
%length(SZ_c_rd)
SZ_c_sum = sum(SZ_c);
SZ_c_diff = Nc - sum(SZ_c);
%return

if Nc - sum(SZ_c) > 0
   SZ_c(Nn) = SZ_c(Nn) + 1 ; 
end


if Nc - sum(SZ_c) > 0
    for iii = 1:Nc-sum(SZ_c)
        jjj = 0;
        while (jjj < 100000)
            new_N = floor((length(my_x)) * rand);
            ctrl_rn = rand;
            if (ctrl_rn < SZ(new_N+1))
                SZ_c(new_N+1) = SZ_c(new_N+1) + 1;
                break
            else
                jjj = jjj + 1;
            end
        end
    end    
end

SZ_c_sum = sum(SZ_c);




%plot(my_x, SZ)
%hold on
%bar(my_x, SZ_c)



% compare 1st, 2nd and 3rd moment of discreze SZ_c and real Schulz-Zimm

% 1st moment "real" SZ 
SZ_1 = Nn;
SZ_1_alt = sum(my_x .* (SZ*Nc)/sum(Nc*SZ)); 

% 1st moment discrete SZ here
SZ_c_1 = sum(my_x .* SZ_c/sum(SZ_c));

dev_1 = (SZ_1 - SZ_c_1)/SZ_1 * 100.;
dev_1_alt = (SZ_1_alt - SZ_c_1)/SZ_1_alt * 100.;
%return

% 2nd moment "real" SZ
SZ_2 = Nn*Nn * (K+1)/K;
SZ_2_alt = sum((my_x.^2) .* (SZ*Nc)/sum(Nc*SZ)); 

% 2nd moment discrete SZ
SZ_c_2 = sum((my_x.^2) .* SZ_c/sum(SZ_c));

dev_2 = (SZ_2 - SZ_c_2) / SZ_2 * 100.;
dev_2_alt = (SZ_2_alt - SZ_c_2)/SZ_2_alt * 100.;

% 3rd moment "real" SZ
SZ_3 = Nn*Nn*Nn *((K+1) * (K+2))/(K*K);
SZ_3_alt = sum((my_x.^3) .* (SZ*Nc)/sum(Nc*SZ)); 

% 3rd moment "discrete" SZ
SZ_c_3 = sum( (my_x.^3) .* SZ_c/sum(SZ_c) );

dev_3 = (SZ_3 - SZ_c_3) / SZ_3 * 100.;
dev_3_alt = (SZ_3_alt - SZ_c_3)/SZ_3_alt * 100.;

% 4th moment "real" SZ
SZ_4 = Nn*Nn*Nn*Nn * ((K+1)*(K+2)*(K+3))/(K*K*K);
SZ_4_alt = sum((my_x.^4) .* (SZ*Nc)/sum(Nc*SZ)); 

% 4th moment discrete SZ
SZ_c_4 = sum( (my_x.^4) .* SZ_c/sum(SZ_c) );

dev_4 = (SZ_4 - SZ_c_4) / SZ_4 * 100.;
dev_4_alt = (SZ_4_alt - SZ_c_4)/SZ_4_alt * 100.;


%if (abs(dev_1) <= 1.0) & (abs(dev_2) <= 1.0) & (abs(dev_3) <= 1.0) & (abs(dev_4) <= 1.0);
%if (abs(dev_1) <= 1.0) & (abs(dev_2) <= 1.0) & (abs(dev_3) <= 1.0) & (sum(SZ_c) == Nc);
%if (abs(dev_1_alt) <= 1.0) & (abs(dev_2_alt) <= 1.0) & (abs(dev_3_alt) <= 1.0) & (abs(dev_4_alt) <= 1.0) & (sum(SZ_c) == Nc) ;
if (abs(dev_1_alt) <= 1.0) & (abs(dev_2_alt) <= 1.0) & (abs(dev_3_alt) <= 1.0) & (abs(dev_4_alt) <= 1.0) & (sum(SZ_c) == Nc) &  (abs(dev_1) <= 1.0) & (abs(dev_2) <= 1.0) & (abs(dev_3) <= 1.0) & (abs(dev_4) <= 1.0);
%if (abs(dev_1_alt) <= 1.0) & (abs(dev_2_alt) <= 1.0) & (abs(dev_3_alt) <= 1.0) & (sum(SZ_c) == Nc) ;
    break
else
    lll=lll+1;
end

end

sum(SZ*Nc)
sum(SZ_c)

% weight averaged chain length
Nw = Nn * (K+1)/K
Nw_c = sum( (my_x.^2) .* SZ_c/sum(SZ_c)) / sum( my_x .* SZ_c/sum(SZ_c))

Nw_over_Nn = Nw / Nn
Nw_c_over_Nn_c = Nw_c / SZ_c_1

figure(1);
clf;
bar(my_x, SZ_c)
hold on
plot(my_x, SZ*Nc, 'r-o')


%% print the distribution out as well into a file

distr_fn= strcat('Schulz_Zimm_Nc',num2str(Nc),'_Navg', num2str(Nmon1+Nmon2), '_sig', num2str(sigma), '_K', num2str(K), '.txt')

fid3 = fopen(distr_fn,'wt');
%fprintf(fid2,'LAMMPS Description of a Polymer Brush with branched Chains\n');
%fprintf(fid2,'branched brush: rho %1.1f final wall distance %3.2f compress for %7.1f steps at t_lj=%2.3f and velocity per wall=%2.3f\n',rho,Lz_final,compress_steps,tlj,vel_per_wall);
fprintf(fid3,'#Polydisperse brush: N_chains = %i, N_avg = %i, sigma = %3.5f K = %3.5f, 1 + 1/K = %3.5f\n',Nc,Nn, sigma,K,polyd);
fprintf(fid3,'#N SZ_sampled SZ_ideal SZ_ideal*Nchains\n');
for blabb = 1:length(my_x)
    fprintf(fid3,'%i %i %3.5f %3.5f\n',my_x(blabb), SZ_c(blabb), SZ(blabb), SZ(blabb)*Nc);
end
fclose(fid3);
%fclose('all');
clear fid3;


%return

%SZ_c_diff = sum(SZ_c) - Nc % some difference is there, either too many or to few




%return


gen=0
if (Nmon2==0)
    gen=1;
elseif ((Nmon2 ~= 0) && (Nmon3 == 0))
    gen=2;
elseif (Nmon3 ~= 0)
    gen=3;
end;

filename = strcat('Poly_Ab_', num2str(Nmon1+Nmon2), '_', num2str(sigma), '_K', num2str(K))
count1=0;
count2=0;
Nangles = 0;
Ndih = 0;
Nimp = 0;





%%%% THIS CALCULATES THE AMOUNT OF SOLVENT PARTICLES FOR A DENSITY OF 0.8

if (solventtype == 2)
    N_solv=floor((rho*Lx*Ly*(Lz_final+2*2^(1/6))-Nmontot*Nc*Nbrush)/2)+1;
elseif (solventtype == 1)   
    N_solv=floor((rho*Lx*Ly*(Lz_final+2*2^(1/6))-Nmontot*Nc*Nbrush)/2)+1;
    %N_solv=2*N_solv;
elseif (solventtype == 0)       
    N_solv=0
end
    
if (Nbrush==1)
    N_solv=0
end


if (Nbrush==1)
    Tatoms = 8;
elseif (Nbrush==2)
    Tatoms = 17;
    if (N_solv == 0)
        Tatoms = Tatoms-1 ;
    end;
    %if (gen == 1)
    %    Tatoms = Tatoms-2 ;
end;

Tbonds = 1;
Tangles = 0;
Tdih = 0;
Timp = 0;


Lx=Lx*sqrt(factor);
Ly=Ly*sqrt(factor);


%Nmontot = sum(my_x .* SZ_c)

%return

%% from chain length distribution generate array of chain lengths and
%% then randomize them
Nmonchain_all = []

for my_idx = [1:1:length(my_x)]
    if SZ_c(my_idx) > 0
       for hhh = [1:1:SZ_c(my_idx)]
           Nmonchain_all = [Nmonchain_all, my_x(my_idx)];
       end
    end
end

%plot([1:length(Nmonchain_all)], Nmonchain_all, 'r-o')

% shuffle randomly
Nmonchain_all = Nmonchain_all(randperm(length(Nmonchain_all)));
%Nmonchain_all = Nmonchain_all(randperm(length(Nmonchain_all)));
%Nmonchain_all = Nmonchain_all(randperm(length(Nmonchain_all)));
%Nmonchain_all(randperm(length(Nmonchain_all)));
%Nmonchain_all(randperm(length(Nmonchain_all)));
%Nmonchain_all(randperm(length(Nmonchain_all)));

%% next aim:
%% swap one element with N = Nn with another chain located at the middle

%% 1st task: find an element with N = Nn
Nn_elements = find(Nmonchain_all == Nn) % indices of these elements
%Nmonchain_all(Nn_elements) % check that they all have N = Nn
my_Nn_element = Nn_elements(randperm(length(Nn_elements),1)) % pick one at random


%% 2nd task: find the middle position on the grid + the element there
%% best is to just guess a few times and then compare later
my_middle_element = floor(Nc/2) + 1%323

%% swap the two elements = x([i j]) = x([j i])
Nmonchain_all([my_Nn_element my_middle_element]) = Nmonchain_all([my_middle_element my_Nn_element])


%hold on
%plot([1:length(Nmonchain_all)], Nmonchain_all, 'b-o')

%return 

%ztot=zeros(Nbrush*Nmontot*Nc,3);%zeros(Nmon1+2*Nmon2,3);
ztot=zeros(Nbrush * sum(my_x .* SZ_c),3) 
type=zeros(90000,1);

for NB=1:Nbrush;
    
    ugrid=zeros(floor(Ly+1),floor(Lx+1));
    lgrid=zeros(floor(Ly+1),floor(Lx+1));
    size(lgrid)
    Ngrid=floor(Ly+1)*floor(Lx+1);
    sc=2^(1/6);
    mz=min(b.z)+sc;
    z=zeros(Nmontot,3);%zeros(Nmon1+2*Nmon2,3);
    kk1=0;kk2=0;
    v=z;
    nn=zeros(90000,3); %%% information on jumping over boundary conditions
    %lower brush coordinates
    
    
    Lx
    Ly
    area=Lx*Ly
    ratio=Lx/Ly
    sqrtNc=sqrt(Nc)
    xgrid=round(sqrtNc*ratio)+4 %-3 %floor(sqrtNc*ratio)
    ygrid=round(sqrtNc/ratio)-3 %+5 %floor(sqrtNc/ratio)
    Nc
    myNc = xgrid*ygrid
    xinc=Lx/xgrid;
    yinc=Ly/ygrid;
    %xgridfloor((myNc-Nc)/2)
    %return
    gridcounter=0
    Nmon1_old = Nmon1
    Nmon2_old = Nmon2
    Nsofar = 0

        for nchains=1:Nc %Nc
            ready = false;
            count=0;
            while ~ready
                count=count+1;
                if (count > 10000)
                    fprintf('Cannot fit polymers on grid. Choose a smaller Nc.\n');
                    return;
                end
                nchains;
                
                jj = (nchains-1)+1-gridcounter*xgrid;
                ii = gridcounter+1;
                if (nchains == xgrid*(gridcounter+1))
                    gridcounter = gridcounter+1;
                end    
                ready=1;
                %return

            end  
            
            
            %get the current chain length
            Nmontot = Nmonchain_all(nchains);
            %Nsofar = 0
            
            %get new Nmon1 and Nmon2
            if Nmontot == 0
                Nmon1 = 0;
                Nmon2 = 0;
                %kk1=kk1+1;
                %graftingpoints(kk1,2) = ii;
                %graftingpoints(kk1,1) = jj;
                %graftingcoords(kk1,1) = (jj-1)*xinc-Lx/2+1+xinc/4; %z(Nsofar+1,2);
                %graftingcoords(kk1,2) = (ii-1)*yinc-Ly/2+1+yinc/4+mod(jj,2)*yinc/2;%z(Nsofar+1,1);                
                "Zero"
                continue
            
            %return
            elseif Nmontot == 1
                Nmon1 = 1;
                Nmon2 = 0;
            elseif Nmontot >= 2
                Nmon1 = 2;
                Nmon2 = Nmontot-2;
            end

            
            
            kk1=kk1+1;
            %z((kk1-1)*Nmontot+1:Nmon1+(kk1-1)*Nmontot,2)=(ii-1)*yinc-Ly/2+1+yinc/4+mod(jj,2)*yinc/2; %; %y coord of polymer
            %z((kk1-1)*Nmontot+1:Nmon1+(kk1-1)*Nmontot,1)=(jj-1)*xinc-Lx/2+1+xinc/4; %+mod(ii,2)*xinc/2; %x coord of polymer
            %z((kk1-1)*Nmontot+1:Nmon1+(kk1-1)*Nmontot,3)=mz-1+[1:Nmon1];  % z coord of polymer
            z(Nsofar+1:Nmon1+Nsofar,2)=(ii-1)*yinc-Ly/2+1+yinc/4+mod(jj,2)*yinc/2; %; %y coord of polymer
            z(Nsofar+1:Nmon1+Nsofar,1)=(jj-1)*xinc-Lx/2+1+xinc/4; %+mod(ii,2)*xinc/2; %x coord of polymer
            z(Nsofar+1:Nmon1+Nsofar,3)=mz-1+[1:Nmon1];  % z coord of polymer
            
            
            graftingpoints(kk1,2) = ii;
            graftingpoints(kk1,1) = jj;
            
            if (Nmon2 > 0) % erste Ordnung Verzweigung
                
                        for yy=1:q
                            yy;
                            z(Nsofar+Nmon1+1+(yy-1)*Nmon2:Nsofar+Nmon1+yy*Nmon2,1)=(jj-1)*xinc-Lx/2+1+xinc/4; %+mod(ii,2)*xinc/2; %+ 0.9*cos((yy-1)/q*2*pi);  %jj-Lx/2+0.9*cos((yy-1)/q*2*pi);                
                            z(Nsofar+Nmon1+1+(yy-1)*Nmon2:Nsofar+Nmon1+yy*Nmon2,2)=(ii-1)*yinc-Ly/2+1+yinc/4+mod(jj,2)*yinc/2; %+0.9*sin((yy-1)/q*2*pi);  %ii-Ly/2+0.9*sin((yy-1)/q*2*pi);%+rl;
                            z(Nsofar+Nmon1+1+(yy-1)*Nmon2:Nsofar+Nmon1+yy*Nmon2,3)=mz+Nmon1-1+[1:Nmon2]; %mz+Nmon1-1.5+[1:Nmon2];

                        end
                        

            end
            
            
            
            
         
%type(1:Nmontot*Nc*NB)=1;            

%middle
type(Nsofar + 1:Nsofar +Nmontot) = 2 + mod(NB+1,2)*6;

%  end
type(Nsofar + Nmontot)  = 3 + mod(NB+1,2)*6;


% head
%type(Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+1,1)=1 + mod(NB+1,2)*6;             
type(Nsofar + 1) = 1 + mod(NB+1,2)*6;

% upper bulk   
%type(Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1+1:Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1+q*Nmon2,1)= 4+mod(NB+1,2)*6; 
%type(Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1+1:Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1+q*Nmon2,1)= 2+mod(NB+1,2)*6;

% lower bulk
%type(Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+2:Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1,1)=2 + mod(NB+1,2)*6;
%type(Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+2:Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1,1)=2 + mod(NB+1,2)*6;

% middle monomer
%type(Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1,1)=2 + mod(NB+1,2)*6; 




% end
%type(Nmontot*Nc*(NB-1)+kk1*Nmontot,1)=5 + mod(NB+1,2)*6;            
%type(Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1+Nmon2:Nmon2:Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1+Nmon2*q,1)=3 + mod(NB+1,2)*6; 

        
            %z((kk1-1)*Nmontot+1:kk1*Nmontot,1)=z((kk1-1)*Nmontot+1:kk1*Nmontot,1)-1;
            %z((kk1-1)*Nmontot+1:kk1*Nmontot,2)=z((kk1-1)*Nmontot+1:kk1*Nmontot,2)-1;
           % z((kk1-1)*Nmontot+1:kk1*Nmontot,1)=z((kk1-1)*Nmontot+1:kk1*Nmontot,1)-1;
           % z((kk1-1)*Nmontot+1:kk1*Nmontot,2)=z((kk1-1)*Nmontot+1:kk1*Nmontot,2)-1;
            z(Nsofar+1:Nsofar+Nmontot,1)=z(Nsofar+1:Nsofar+Nmontot,1)-1;
            z(Nsofar+1:Nsofar+Nmontot,2)=z(Nsofar+1:Nsofar+Nmontot,2)-1;
            
            graftingcoords(kk1,1) = z(Nsofar+1,2);
            graftingcoords(kk1,2) = z(Nsofar+1,1);

               %type((kk1)*Nmontot-1,1)=4; 
               
            Nsofar = Nsofar + Nmontot;
        end
        %return
        
        % second brush on top
        %ztot(Nmontot*Nc*(NB-1)+1:NB*Nmontot*Nc,1:2)=z(1:Nmontot*Nc,1:2);
        %ztot(Nmontot*Nc*(NB-1)+1:NB*Nmontot*Nc,3)=z(1:Nmontot*Nc,3)*(-1)^(mod(NB+1,2))+ mod(Nbrush+1,2)*Lz/2*(-1)^(mod(NB,2));
        ztot(Nsofar*(NB-1)+1:Nsofar*NB, 1:2) = z(1:Nsofar,1:2);
        ztot(Nsofar*(NB-1)+1:Nsofar*NB, 3) = z(1:Nsofar,3)*(-1)^(mod(NB+1,2))+ mod(Nbrush+1,2)*Lz/2*(-1)^(mod(NB,2));
end


%return
%%% GENERATE SOLVENT MOLECULES



y=zeros(2,3);
xjump=0;
yjump=0;
zjump=0;
nlayers=0;
for i=1:N_solv
        y((i-1)*2+1:i*2,1) = -Lx/2 + [1:2] + xjump +0.2;
        y((i-1)*2+1:i*2,2) = -Ly/2 + yjump+0.2;
        y((i-1)*2+1:i*2,3) = zjump;
        xjump=xjump+3;
        if (xjump > Lx-1)
            yjump=yjump+2;
            xjump=0;
        end
        if (yjump > Ly-1)
            yjump=0;
            nlayers=nlayers+1;
            zjump=zjump + 3*nlayers*(-1)^(nlayers+1);
        end
%         if (zjump == 6)
%             yjump=0
%             zjump=zjump-9;
%         end        
end;


Lx=Lx/sqrt(factor);  %%% back to original size
Ly=Ly/sqrt(factor);  %%% back to original size


typeA= Nbrush * sum(my_x .* SZ_c); %Nbrush*Nc*(Nmontot);
typeS = 2*N_solv;

type(typeA+1:typeA+typeS)=17;  % type of solvent

%typeA=(kk1+kk2)*(Nmontot);    % # of monomers in the brush (includes grafted heads)
typeA1=Nbrush*(Nmontot-1)*Nc;
typeB1=Nbrush;
%typeS = 2*N_solv;


%x=zeros(size(z,1)+typeS+number1/2,3);
x=zeros(size(ztot,1)+typeS+nnn*Nbrush,3);


%number3=size(z,1)+number1/2+typeS; %size(z,1)+number1+typeS;
number3=size(ztot,1)+nnn*Nbrush+typeS;

xcorrection=floor(Lx*sqrt(factor)+1)/Lx; %floor(Lx*sqrt(factor)+1)/sqrt(factor)/Lx
ycorrection=floor(Ly*sqrt(factor)+1)/Ly; %floor(Ly*sqrt(factor)+1)/sqrt(factor)/Ly

xcorrection=1.0;
ycorrection=1.0;

% z(:,1)=z(:,1)./sqrt(factor)/xcorrection; %%% rescale x and y coordinates in case of stretching system
% z(:,2)=z(:,2)./sqrt(factor)/ycorrection;
ztot(:,1)=ztot(:,1)./xcorrection; %%% rescale x and y coordinates in case of stretching system
ztot(:,2)=ztot(:,2)./ycorrection;

x(1:typeA,1:3)=ztot;

% x(1:typeA,1:3)
% size(x)
% number3-typeA
% number1/2

if (N_solv > 0)
    x(typeA+1:typeA+typeS,1:3)=y;
end

% size(z)
% size(a.x)
% size(x)

x(typeA+typeS+1:typeA+typeS+nnn,1)=b.x;
x(typeA+typeS+1:typeA+typeS+nnn,2)=b.y;
x(typeA+typeS+1:typeA+typeS+nnn,3)=b.z-Lz/2*mod(Nbrush+1,2);

if (Nbrush == 2)
    x(typeA+typeS+nnn+1:typeA+typeS+Nbrush*nnn,1)=b.x;
    x(typeA+typeS+nnn+1:typeA+typeS+Nbrush*nnn,2)=b.y;
    x(typeA+typeS+nnn+1:typeA+typeS+Nbrush*nnn,3)=b.z+Lz/2;
end


type(typeA+typeS+1:typeA+typeS+nnn) = 4;    % lower wall
if (Nbrush == 2)
    type(typeA+typeS+nnn+1:typeA+typeS+Nbrush*nnn) = 12;    % upper wall
end
%type(typeA+typeS+number1/2+1:typeA+typeS+number1) = 8;    % upper wall

vel_per_wall=0.01;
tlj=0.0025;
Lzbegin=Lz;
compress_steps = (Lzbegin-(Lz_final+2*2^(1/6)))/(vel_per_wall*tlj*2);

compress_steps1 = (Lzbegin-(Lz_final1+2*2^(1/6)))/(vel_per_wall*tlj*2);
compress_steps2 = (Lzbegin-(Lz_final2+2*2^(1/6)))/(vel_per_wall*tlj*2);
compress_steps3 = (Lzbegin-(Lz_final3+2*2^(1/6)))/(vel_per_wall*tlj*2);


lmp_filename=strcat(filename, 'huge.data');

fid2 = fopen(lmp_filename,'wt');
%fprintf(fid2,'LAMMPS Description of a Polymer Brush with branched Chains\n');
fprintf(fid2,'branched brush: rho %1.1f final wall distance %3.2f compress for %7.1f steps at t_lj=%2.3f and velocity per wall=%2.3f\n',rho,Lz_final,compress_steps,tlj,vel_per_wall);
fprintf(fid2,'\n');

Natoms = number3;
N_solv;

%% update this here, too

if (solventtype == 2)
    %Nbonds = Nbrush*Nc*(Nmontot-1)+N_solv       %% f�r 2-fach Verzweigung ist es Nmontot-1    %(kk1+kk2)*(Nmontot-1)+N_solv
    %Nbonds = Nbrush * sum((my_x -1) .* SZ_c) + Nsolv;
    Nbonds = Nbrush * sum((my_x(2:end) -1) .* SZ_c(2:end)) + Nsolv;
else
    %Nbonds = Nbrush * sum((my_x -1) .* SZ_c);     %Nbrush*Nc*(Nmontot-1)
    Nbonds = Nbrush * sum((my_x(2:end) -1) .* SZ_c(2:end));
end

fprintf(fid2,'%d atoms\n',Natoms);
fprintf(fid2,'%d bonds\n',Nbonds);
fprintf(fid2,'%d angles\n',Nangles);
fprintf(fid2,'%d dihedrals\n',Ndih);
fprintf(fid2,'%d impropers\n',Nimp);
fprintf(fid2,'\n');
       
fprintf(fid2,'%d atom types\n',Tatoms);
fprintf(fid2,'%d bond types\n',Tbonds);
fprintf(fid2,'%d angle types\n',Tangles);
fprintf(fid2,'%d dihedral types\n',Tdih);
fprintf(fid2,'%d improper types\n',Timp);
fprintf(fid2,'\n');

%Box measurements
%for periodic systems this is box size
%for non-periodic it is min/max extent of atoms
fprintf(fid2,'%6.5f %6.5f xlo xhi\n',-Lx/2, Lx/2);
fprintf(fid2,'%6.5f %6.5f ylo yhi\n',-Ly/2, Ly/2);
%fprintf(fid2,'%6.5f %6.5f zlo zhi\n',-Lz/2, Lz/2);
%fprintf(fid2,'%6.5f %6.5f zlo zhi\n',0.00-mod(Nbrush+1,2)*Lz/2-3, Lz-mod(Nbrush+1,2)*Lz/2+3);
fprintf(fid2,'%6.5f %6.5f zlo zhi\n',0.00-mod(Nbrush+1,2)*Lz/2-3, my_x(find(SZ_c,1,'last'))-mod(Nbrush+1,2)*Lz/2+3);
fprintf(fid2,'\n');
%Masses

 % 1 mass
 % ...
 % N mass                           (N = # of atom types)

 %Masses - same number as atom types
fprintf(fid2,'Masses\n');
fprintf(fid2,'\n'); 
fprintf(fid2,'%d %1.2f\n',1, 1.0);
fprintf(fid2,'%d %1.2f\n',2, 1.0);
fprintf(fid2,'%d %1.2f\n',3, 1.0);
fprintf(fid2,'%d %1.2f\n',4, 1.0);
fprintf(fid2,'%d %1.2f\n',5, 1.0);
fprintf(fid2,'%d %1.2f\n',6, 1.0);
if (Tatoms >= 7)
fprintf(fid2,'%d %1.2f\n',7, 1.0);
end;
if (Tatoms >= 8)
fprintf(fid2,'%d %1.2f\n',8, 1.0);
end;
if (Tatoms >= 9)
    fprintf(fid2,'%d %1.2f\n',9, 1.0);
end;
if (Tatoms >= 10)
    fprintf(fid2,'%d %1.2f\n',10, 1.0);
end;
if (Tatoms >= 11)
    fprintf(fid2,'%d %1.2f\n',11, 1.0);
end;
if (Tatoms >= 12)
    fprintf(fid2,'%d %1.2f\n',12, 1.0);
end;
if (Tatoms >= 13)
    fprintf(fid2,'%d %1.2f\n',13, 1.0);
end;
if (Tatoms >= 14)
    fprintf(fid2,'%d %1.2f\n',14, 1.0);
end;
if (Tatoms >= 15)
    fprintf(fid2,'%d %1.2f\n',15, 1.0);
    
end;
if (Tatoms >= 16)
    fprintf(fid2,'%d %1.2f\n',16, 1.0);
end;
if (Tatoms >= 17)
    fprintf(fid2,'%d %1.2f\n',17, 1.0);
end;
fprintf(fid2,'\n');
 

%Atoms
% n molecule-tag atom-type q x y z nx ny nz
fprintf(fid2,'Atoms\n');
fprintf(fid2,'\n'); 
%molec_dummy = 1;
molec_dummy = 1; % counts consecutive, jumps over chains of length 0
molec_counter = 0; % counts til Nc
%atom_dummy = 1;
type_dummy = 1;
%number3
%typeB1*Nmon
%typeA/Nmon*(Nmon-1)


% save resids of molecules with length 0 and 1 to check if this stuff works
ids_0 = [];
ids_1 = [];

Nsofar = 0;


Nmontot = 0; %Nmonchain_all(molec_dummy);
while Nmontot == 0
    %molec_dummy = molec_dummy;
    molec_counter = molec_counter +1;
    Nmontot = Nmonchain_all(molec_counter);
    if Nmontot == 0
        "ZERO"
        ids_0 = [ids_0, molec_dummy]
    end
end

for i=1:number3
    %fprintf(fid2,'%i %i %i %6.5f %6.5f %6.5f %i %i %i\n',i, molec_dummy, type(i), x(i,1),x(i,2),x(i,3),nn(i,1),nn(i,2),nn(i,3));
    fprintf(fid2,'%i %i %i %6.5f %6.5f %6.5f\n',i, molec_dummy, type(i), x(i,1),x(i,2),x(i,3));
    if (Nmontot == 0)
        molec_dummy = molec_dummy -1;
    end
    
    %if ((mod(i,Nmontot) == 0) && (i < typeA+1))
    if ((i >= Nsofar+Nmontot) && (i < typeA+1))    
        molec_dummy = molec_dummy + 1;
        molec_counter = molec_counter +1;
        Nsofar = Nsofar + Nmontot;
        if molec_counter <= length(Nmonchain_all)
            Nmontot = Nmonchain_all(molec_counter);
            if Nmontot == 0
                "zzzero"
                ids_0 = [ids_0, molec_dummy]
                %return
            elseif Nmontot == 1
                "ooone"
                ids_1 = [ids_1, molec_dummy]
            end
        end
    end
    
    if (solventtype == 2)
        if ((i > typeA) && (i < typeA+typeS+1) && (mod(i-Nmontot,2)==0))
            molec_dummy = molec_dummy + 1;
        end
    elseif (solventtype == 1)
        if ((i > typeA) && (i < typeA+typeS+1))
            molec_dummy = molec_dummy + 1;
        end     
    end
    
    if (i > typeA+typeS) 
        molec_dummy = molec_dummy + 1;
    end
end
fprintf(fid2,'\n');

% Velocities if needed

% Bonds
% n bond-type atom-1 atom-2
fprintf(fid2,'Bonds\n');
fprintf(fid2,'\n');
bond_index = 1;
bond_type = 1;
%for i=1:typeA/Nmon1
%      for j=1:Nmon1-1
%         fprintf(fid2,'%i %i %i %i\n',bond_index, bond_type, (i-1)*Nmon1+j,(i-1)*Nmon1+j+1);
%         bond_index = bond_index + 1;
%      end
%end




Nsofar = 0;
for i=1:Nc*Nbrush% loop over polymer chains
    %for j=1:Nmontot-1 % loop over
    Nmontot=Nmonchain_all(i);
    
     if Nmontot == 0
         Nmon1 = 0;
         Nmon2 = 0;
         continue
     elseif Nmontot == 1
         Nmon1 = 1;
         Nmon2 = 0;
         Nsofar = Nsofar +1;
         continue
     elseif Nmontot >= 2
         Nmon1 = 2;
         Nmon2 = Nmontot-2;
     end
    
    for j=1:Nmon1-1
        %1stgen and middle arm
        %if (j < Nmon1)
            %fprintf(fid2,'%i %i %i %i\n', bond_index, bond_type, (i-1)*Nmontot+j,(i-1)*Nmontot+j+1);
            fprintf(fid2,'%i %i %i %i\n', bond_index, bond_type, Nsofar+j, Nsofar+j+1);
            bond_index = bond_index + 1;
        %end
    end
    for yy=1:q
        for j=1:Nmon2
        %for yy=1:q
            %2ndgen arms
            if (j == 1) %Nmon1+Nmon2*()
                fprintf(fid2,'%i %i %i %i\n', bond_index, bond_type, Nsofar+Nmon1,Nsofar+Nmon1+(yy-1)*Nmon2+j);
                bond_index = bond_index + 1;
            else
                fprintf(fid2,'%i %i %i %i\n', bond_index, bond_type, Nsofar+Nmon1+(yy-1)*Nmon2+j-1,Nsofar+Nmon1+(yy-1)*Nmon2+j);
                bond_index = bond_index + 1;
            end
            
            %end        
            %if ((j > Nmon1+Nmon2*(yy-1)) && (j < Nmon1+2*Nmon2))
            %    fprintf(fid2,'%i %i %i %i\n', bond_index, bond_type, (i-1)*Nmontot+j,(i-1)*Nmontot+j+1);
            %    bond_index = bond_index + 1;
            %end
        end 
    end
    
    Nsofar = Nsofar + Nmontot;
end    

if (solventtype == 2)
    for i=1:N_solv
          fprintf(fid2,'%i %i %i %i\n',bond_index, bond_type, typeA+(i-1)*2+1,typeA+(i-1)*2+2);
          bond_index = bond_index + 1;
    end
end

kk1;
kk2;

fclose(fid2);
%fclose('all');
clear fid2;


toc
figure(10)
clf;
plot3(x(:,1),x(:,2),x(:,3),'.')
% min(graftingcoords(:,1))
% min(graftingcoords(:,2))
% max(graftingcoords(:,1))
% max(graftingcoords(:,2))

density_after_compression = (Nc*Nmontot*2+N_solv*2)/(Lx*Ly*(Lz_final+2*2^(1/6)));
Lx;
Ly;
border=(Lz_final+2*2^(1/6))/2;
Lz_final;
Lz_final=Lz_final+2*2^(1/6);
v=0.01;
dt=0.0025;
csteps=(Lz-Lz_final)/(2*v*dt);


SIGMA_FINAL=Nc/(Lx*Ly)

SZ_c(1:2)
longest_chain = my_x(find(SZ_c,1,'last'))
