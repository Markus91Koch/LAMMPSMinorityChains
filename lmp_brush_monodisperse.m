
% make a read file for LAMMPS
clear;
tic
%numbers

Nmon1 = 2 %36;
Nmon2 = 128 - Nmon1 %35;
Nmon3 = 0;

q = 1;

Lz = 290;

Nbrush = 1; %% TWO BRUSHES OR 1

solventtype = 0; %% no solvent (0), monomeric (1) or dimeric solvent (2)

if ((solventtype ~=1) && (solventtype ~= 0) && (solventtype ~= 2)) 
    solventtype
    fprintf('Choose value 0 (none), 1 (monomeric) or 2 (dimeric) for solventtype!\n'); 
    return;
end   

%%% gut: "factor" ist quadrat einer natuerlichen zahl oder
%%% eines bruchs natürlicher zahlen
factor=1; %4   %16/9; %25/9; %4; %% increase initial density by this factor

Nmontot = Nmon1+Nmon2*q %+Nmon3*4

sigma_final = 0.100; %0.025 %2

sigma = sigma_final / factor;

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

rho = 0.85;
Lz_final = 12.0; %100.00  %12.0%17.5;

Lz_final1 = 12.0;
Lz_final2 = 14.25;
Lz_final3 = 17.5;

xrows=(152*3/4-14)*1.5;
yrows=(32*3/4+6)*1.5;

nnn = xrows*yrows; %% number of lattice atoms

%% array to save coordinates of l21;attice 
b = struct('x',zeros(nnn,1),'y',zeros(nnn,1),'z',zeros(nnn,1));

%% increment of lattice 
dx=0.52500000;
dy=0.9093250000;

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
Lx=xmax-xmin+dx
Ly=ymax-ymin+dy
Nc = floor(sigma*Lx*Ly*factor)+1

 Nc = Nc -1 + 0
%Nc = floor(sigma*Lx_final*Ly_final*factor)+1
%Nc=2

% alllocate space for grafting coords etc.
graftingpoints = zeros(2*Nc,2);
graftingcoords = zeros(2*Nc,2);


gen=0
if (Nmon2==0)
    gen=1;
elseif ((Nmon2 ~= 0) && (Nmon3 == 0))
    gen=2;
elseif (Nmon3 ~= 0)
    gen=3;
end;


%% devise filename
filename = strcat('Ab_', num2str(Nmon1+Nmon2), '_', num2str(sigma))

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
end;

Tbonds = 1;
Tangles = 0;
Tdih = 0;
Timp = 0;


Lx = Lx*sqrt(factor);
Ly = Ly*sqrt(factor);


ztot = zeros(Nbrush*Nmontot*Nc,3);%zeros(Nmon1+2*Nmon2,3);
type = zeros(90000,1);

for NB=1:Nbrush;
    
    ugrid = zeros(floor(Ly+1),floor(Lx+1));
    lgrid = zeros(floor(Ly+1),floor(Lx+1));
    size(lgrid)
    Ngrid = floor(Ly+1)*floor(Lx+1);
    sc = 2^(1/6);
    mz = min(b.z)+sc;
    z = zeros(Nmontot,3);%zeros(Nmon1+2*Nmon2,3);
    kk1 = 0;
    kk2 = 0;
    v = z;
    nn = zeros(90000,3); %%% information on jumping over boundary conditions
    %lower brush coordinates
    
    Lx
    Ly
    area = Lx*Ly
    ratio = Lx/Ly
    sqrtNc = sqrt(Nc)
    xgrid = round(sqrtNc*ratio)+4 %-3 %floor(sqrtNc*ratio)
    ygrid = round(sqrtNc/ratio)-3 %+5 %floor(sqrtNc/ratio)
    Nc
    myNc = xgrid*ygrid
    xinc = Lx/xgrid;
    yinc = Ly/ygrid;

    gridcounter = 0

        for nchains=1:Nc %Nc
            ready = false;
            count = 0;
            while ~ready
                count=count+1;
                
                %%% break condition of loop
                if (count > 10000)
                    fprintf('Cannot fit polymers on grid. Choose a smaller Nc.\n');
                    return;
                end

                nchains
                
                jj = (nchains-1)+1-gridcounter*xgrid
                ii = gridcounter+1
                if (nchains == xgrid*(gridcounter+1))
                    gridcounter = gridcounter+1;
                end    
                ready=1;

            end  
            
            
            kk1=kk1+1;
            z((kk1-1)*Nmontot+1:Nmon1+(kk1-1)*Nmontot,2) = (ii-1)*yinc - Ly/2 + 1 + yinc/4 + mod(jj,2)*yinc/2; %; %y coord of polymer
            z((kk1-1)*Nmontot+1:Nmon1+(kk1-1)*Nmontot,1) = (jj-1)*xinc - Lx/2 + 1 + xinc/4; %+mod(ii,2)*xinc/2; %x coord of polymer
            z((kk1-1)*Nmontot+1:Nmon1+(kk1-1)*Nmontot,3) = mz - 1 + [1:Nmon1];  % z coord of polymer
            
            graftingpoints(kk1,2) = ii;
            graftingpoints(kk1,1) = jj;
            
            if (Nmon2 > 0) % erste Ordnung Verzweigung
                
                        for yy=1:q
                            yy;
                            z((kk1-1)*Nmontot+Nmon1+1+(yy-1)*Nmon2:(kk1-1)*Nmontot+Nmon1+yy*Nmon2,1) = (jj-1) * xinc - Lx/2 + 1 + xinc/4;
                            z((kk1-1)*Nmontot+Nmon1+1+(yy-1)*Nmon2:(kk1-1)*Nmontot+Nmon1+yy*Nmon2,2) = (ii-1) * yinc - Ly/2 + 1 + yinc/4 + mod(jj,2) * yinc/2;
                            z((kk1-1)*Nmontot+Nmon1+1+(yy-1)*Nmon2:(kk1-1)*Nmontot+Nmon1+yy*Nmon2,3) = mz + Nmon1 - 1 + [1:Nmon2];

                        end
                        
            end
            
            
            %%%% assign atom types          
            % head
            type(Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+1,1) = 1 + mod(NB+1,2) * 6;
      
            % upper bulk   
            type(Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1+1:Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1+q*Nmon2,1) = 2+mod(NB+1,2) * 6;

            % lower bulk
            type(Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+2:Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1,1) = 2 + mod(NB+1,2) * 6;

            % middle monomer
            type(Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1,1) = 2 + mod(NB+1,2) * 6;

            % end
            type(Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1+Nmon2:Nmon2:Nmontot*Nc*(NB-1)+(kk1-1)*Nmontot+Nmon1+Nmon2*q,1) = 3 + mod(NB+1,2) * 6; 

            z((kk1-1)*Nmontot+1:kk1*Nmontot,1) = z((kk1-1)*Nmontot+1:kk1*Nmontot,1) - 1;
            z((kk1-1)*Nmontot+1:kk1*Nmontot,2) = z((kk1-1)*Nmontot+1:kk1*Nmontot,2) - 1;
            
            graftingcoords(kk1,1) = z((kk1-1)*Nmontot+1,2);
            graftingcoords(kk1,2) = z((kk1-1)*Nmontot+1,1);
            
        end

        ztot(Nmontot*Nc*(NB-1)+1:NB*Nmontot*Nc,1:2) = z(1:Nmontot*Nc,1:2);
        ztot(Nmontot*Nc*(NB-1)+1:NB*Nmontot*Nc,3) = z(1:Nmontot*Nc,3)*(-1)^(mod(NB+1,2))+ mod(Nbrush+1,2)*Lz/2*(-1)^(mod(NB,2));

end

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
            yjump = yjump+2;
            xjump = 0;
        end
        if (yjump > Ly-1)
            yjump=0;
            nlayers = nlayers+1;
            zjump = zjump + 3*nlayers*(-1)^(nlayers+1);
        end
end;


Lx=Lx/sqrt(factor);  %%% back to original size
Ly=Ly/sqrt(factor);  %%% back to original size


typeA=Nbrush*Nc*(Nmontot);
typeS = 2*N_solv;

type(typeA+1:typeA+typeS)=17;  % type of solvent

%typeA=(kk1+kk2)*(Nmontot);    % # of monomers in the brush (includes grafted heads)
typeA1=Nbrush*(Nmontot-1)*Nc;
typeB1=Nbrush;
%typeS = 2*N_solv;


x=zeros(size(ztot,1)+typeS+nnn*Nbrush,3);

number3=size(ztot,1)+nnn*Nbrush+typeS;

xcorrection=floor(Lx*sqrt(factor)+1)/Lx  %floor(Lx*sqrt(factor)+1)/sqrt(factor)/Lx
ycorrection=floor(Ly*sqrt(factor)+1)/Ly  %floor(Ly*sqrt(factor)+1)/sqrt(factor)/Ly

xcorrection=1.0
ycorrection=1.0

ztot(:,1)=ztot(:,1)./xcorrection; %%% rescale x and y coordinates in case of stretching system
ztot(:,2)=ztot(:,2)./ycorrection;

x(1:typeA,1:3)=ztot;

if (N_solv > 0)
    x(typeA+1:typeA+typeS,1:3)=y;
end

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

vel_per_wall = 0.01;
tlj = 0.0025;
Lzbegin = Lz;
compress_steps = (Lzbegin-(Lz_final+2*2^(1/6)))/(vel_per_wall*tlj*2);

compress_steps1 = (Lzbegin-(Lz_final1+2*2^(1/6)))/(vel_per_wall*tlj*2);
compress_steps2 = (Lzbegin-(Lz_final2+2*2^(1/6)))/(vel_per_wall*tlj*2);
compress_steps3 = (Lzbegin-(Lz_final3+2*2^(1/6)))/(vel_per_wall*tlj*2);


%%% Write LAMMPS file

lmp_filename = strcat(filename, 'huge.data')

fid2 = fopen(lmp_filename,'wt');
%fprintf(fid2,'LAMMPS Description of a Polymer Brush with branched Chains\n');
fprintf(fid2,'branched brush: rho %1.1f final wall distance %3.2f compress for %7.1f steps at t_lj=%2.3f and velocity per wall=%2.3f\n',rho,Lz_final,compress_steps,tlj,vel_per_wall);
fprintf(fid2,'\n');

Natoms = number3
N_solv

if (solventtype == 2)
    Nbonds = Nbrush*Nc*(Nmontot-1)+N_solv       %% f�r 2-fach Verzweigung ist es Nmontot-1    %(kk1+kk2)*(Nmontot-1)+N_solv
else
    Nbonds = Nbrush*Nc*(Nmontot-1)
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
fprintf(fid2,'%6.5f %6.5f zlo zhi\n',0.00-mod(Nbrush+1,2)*Lz/2-3, Lz-mod(Nbrush+1,2)*Lz/2+3);
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
molec_dummy = 1;
%atom_dummy = 1;
type_dummy = 1;
%number3
%typeB1*Nmon
%typeA/Nmon*(Nmon-1)
for i=1:number3
    %fprintf(fid2,'%i %i %i %6.5f %6.5f %6.5f %i %i %i\n',i, molec_dummy, type(i), x(i,1),x(i,2),x(i,3),nn(i,1),nn(i,2),nn(i,3));
    fprintf(fid2,'%i %i %i %6.5f %6.5f %6.5f\n',i, molec_dummy, type(i), x(i,1),x(i,2),x(i,3));
    if ((mod(i,Nmontot) == 0) && (i < typeA+1))
        molec_dummy = molec_dummy + 1;
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

% Write Bonds
% n bond-type atom-1 atom-2
fprintf(fid2,'Bonds\n');
fprintf(fid2,'\n');
bond_index = 1;
bond_type = 1;

for i=1:Nc*Nbrush% loop over polymer chains
    %for j=1:Nmontot-1 % loop over
    for j=1:Nmon1-1
        %1stgen and middle arm
        %if (j < Nmon1)
            fprintf(fid2,'%i %i %i %i\n', bond_index, bond_type, (i-1)*Nmontot+j,(i-1)*Nmontot+j+1);
            bond_index = bond_index + 1;
        %end
    end
    for yy=1:q
        for j=1:Nmon2
        %for yy=1:q
            %2ndgen arms
            if (j == 1) %Nmon1+Nmon2*()
                fprintf(fid2,'%i %i %i %i\n', bond_index, bond_type, (i-1)*Nmontot+Nmon1,(i-1)*Nmontot+Nmon1+(yy-1)*Nmon2+j);
                bond_index = bond_index + 1;
            else
                fprintf(fid2,'%i %i %i %i\n', bond_index, bond_type, (i-1)*Nmontot+Nmon1+(yy-1)*Nmon2+j-1,(i-1)*Nmontot+Nmon1+(yy-1)*Nmon2+j);
                bond_index = bond_index + 1;
            end
            
        end 
    end
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
plot3(x(:,1),x(:,2),x(:,3),'.')

density_after_compression = (Nc*Nmontot*2+N_solv*2)/(Lx*Ly*(Lz_final+2*2^(1/6)))
Lx
Ly
border=(Lz_final+2*2^(1/6))/2
Lz_final
Lz_final=Lz_final+2*2^(1/6)
v=0.01;
dt=0.0025;
csteps=(Lz-Lz_final)/(2*v*dt)


SIGMA_FINAL=Nc/(Lx*Ly)
