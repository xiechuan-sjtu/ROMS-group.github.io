%  Interpolate Glosea5 data to ROMS model open boundaries 
%  Ross sea and Amundsen Sea Model (RAS)
%  Written by Xie Chuan
%  Create date: 2021-9-10
%
clear all; 
close all; 
clc;

casename='RAS';
glosea5_data_path='G:\Ross_amundsen_roms_model\boundary_file\B1_Glosea5_before_interp\'; % path of data file 
glosea5_grid_path='G:\Ross_amundsen_roms_model\boundary_file\A1_Glosea5_rawdata\Grid\'; % path of Glosea5 grid data file 
outdata_path='G:\Ross_amundsen_roms_model\boundary_file\B2_Glosea5_after_interp\'; % path of data file 
ROMS_grid_path='G:\Ross_amundsen_roms_model\grid_file\'; % path of ROMS grid 
% ----------------------------- input file --------------------------------
ROMS_grid=[ROMS_grid_path casename '_grd_32layer_new.nc']; % ROMS grid file
Glosea5_grid=[glosea5_grid_path 'Glosea5_grid_ross.mat'];
% ----------------- define time period -----------------
data_year = 1998;  % start year   
nf = 12;
%---------------------------------------------------------------------------------
%
% ---------------------- data file ----------------------
datafile=[];
outfile=[]; % output file name
for ik=1:nf
    if ik < 10
        tk=[num2str(0)  num2str(ik)];
    elseif ik > 9 && ik < 100
        tk=[ num2str(ik)];
    else
        tk=num2str(ik);
    end
    datafile0=[glosea5_data_path 'Glosea5_' num2str(data_year) tk '.mat'];
    datafile=[datafile;datafile0];
    outfile0=[outdata_path casename '_bry_data_' num2str(data_year) tk '.mat']; % output file name
    outfile=[outfile;outfile0];
end

% ------------------ S-coordinate setup ------------------
%
N = 32;           % No. of vertical layers
Vtransform = 2;   % Vertical transformation equation: 1, original; 2, new
Vstretching = 4;  % Vstretching   Vertical stretching function:
%                    Vstretching = 1,  original (Song and Haidvogel, 1994)
%                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS, 2005)
%                    Vstretching = 3,  R. Geyer BBL refinement
%                    Vstretching = 4,  A. Shchepetkin (UCLA-ROMS, 2010)
theta_s = 4;      % S-coordinate surface control parameter (scalar)
theta_b = 2;      % S-coordinate bottom control parameter (scalar)
hc = 200;          % Width (m) of surface or bottom boundary layer in which
                  %  higher vertical resolution is required during
                  %  stretching (scalar)
igrid = 0;        % Depth grid type logical switch:
%                   igrid=1  => density points
%                   igrid=2  => streamfunction points
%                   igrid=3  => u-velocity points
%                   igrid=4  => v-velocity points
%                   igrid=5  => w-velocity points
% 
% --- select open boundaries (by true/false) ---
obc_east=true;
obc_west=true;
obc_north=true;
obc_south=false;
%
% -------------------------------------------------------------------------
% ------------------------------ load data --------------------------------
% -------------------------------------------------------------------------
disp('  Loading grid data ...');
% ---- RHO grid ----
lon_rho=ncread(ROMS_grid,'lon_rho')';
lat_rho=ncread(ROMS_grid,'lat_rho')';
mask_rho=ncread(ROMS_grid,'mask_rho')';

% ---- u grid ----
lon_u=ncread(ROMS_grid,'lon_u')';
lat_u=ncread(ROMS_grid,'lat_u')';
mask_u=ncread(ROMS_grid,'mask_u')';

% ---- v grid ----
lon_v=ncread(ROMS_grid,'lon_v')';
lat_v=ncread(ROMS_grid,'lat_v')';
mask_v=ncread(ROMS_grid,'mask_v')';

% ---- angle ----
ang = ncread(ROMS_grid,'angle')';

% ---- depth ----
h = ncread(ROMS_grid,'h');

% ---- ice shelf ----
zeta_iceshelf = ncread(ROMS_grid,'zice');    % 考虑冰架的影响

% ---- grid dimensions ----
[Mp Lp]=size(lon_rho);   % Number of exterior RHO-points (Lp -> I, Mp -> J);
Lm=Lp-2;                 % Number of I-direction INTERIOR RHO-points
Mm=Mp-2;                 % Number of J-direction INTERIOR RHO-points
L = Lp-1;
M = Mp-1;

% -----------------------------------------consider ice shelf)----------------
% -------------------------------------------------------------------------
disp('  Calculating z level ...');

% for rho points
igrid=1; 
zr=set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N,igrid,h,zeta_iceshelf,0);
zr=permute(zr,[3,2,1]); % rearrange dimensions 
 % for u points
igrid=3;
zu=set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N,igrid,h,zeta_iceshelf,0);
zu=permute(zu,[3,2,1]); % rearrange dimensions 
% for v points
igrid=4; 
zv=set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N,igrid,h,zeta_iceshelf,0);
zv=permute(zv,[3,2,1]); % rearrange dimensions 
% for w points
igrid=5; 
zw=set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N,igrid,h,zeta_iceshelf,0);
zw=permute(zw,[3,2,1]); % rearrange dimensions 
%
h=h';
% -------------------------------------------------------------------------
% ------------------------ Generate boundary data -------------------------
% -------------------------------------------------------------------------
disp(' Gernerating boundary data ...');
%
% ----------- load data ------------
%
load(Glosea5_grid);
lons=lon_Glosea5';
lats=lat_Glosea5';

for fk=1:12
     disp(['    => Processing boundary data file ' num2str(fk) ' of ' num2str(nf)]);
     load(datafile(fk,:));
     salts=salinity; clear salinity;
     temps=temp; clear temp;
     sshs=ssh; clear ssh;
     us=u; clear u;
     vs=v; clear v;
     
     salts=permute(salts,[4 3 2 1]);
     temps=permute(temps,[4 3 2 1]);
     sshs=permute(sshs,[3 2 1]);
     us=permute(us,[4 3 2 1]);
     vs=permute(vs,[4 3 2 1]);
%
     sd=-1*flipud(depth);    % flip depth and change sign
%
     time_date=td; % date of data
     tt=time_date;
%
     [nt nz ny nx]=size(salts);   % size of data
% 
     lon_data=reshape(lons,nx*ny,1);  % change lon and lat to 1-D array
     lat_data=reshape(lats,nx*ny,1);
%
% np=9; % NO. of nearest data points that used to interpolate NaN data
% 
     if obc_east
        zeta_east=[];
        salt_east=[];
        temp_east=[];
        u_east=[];
        v_east=[];
     end

     if obc_west
        zeta_west=[];
        salt_west=[];
        temp_west=[];
        u_west=[];
        v_west=[];
     end

     if obc_north
        zeta_north=[];
        salt_north=[];
        temp_north=[];
        u_north=[];
        v_north=[];
     end

     if obc_south
        zeta_south=[];
        salt_south=[];
        temp_south=[];
        u_south=[];
        v_south=[];
     end
%
% --------- interpolation ---------
     disp('      -> Interpolating boundary ...');

     for tk=1:nt
        disp(['        : Time step ' num2str(tk) ' of ' num2str(nt)]);
    %
% -------- water surface level ---------
        tem0=squeeze(sshs(tk,:,:));
        tem0=reshape(tem0,nx*ny,1);
        tem=[lon_data,lat_data,tem0]; % combine location and data together
        % remove NaN
        kk=find(isnan(sum(tem')));
        tem(kk,:)=[]; 
        tem=double(tem);
        % - interpolaation -
        F = scatteredInterpolant(tem(:,1),tem(:,2),tem(:,3),'linear','nearest');
        % --- eastern boundary ---
        if obc_east
            obc_lon=lon_rho(:,end);
            obc_lat=lat_rho(:,end);
            obc_mask=mask_rho(:,end);
            temi=F(obc_lon,obc_lat); % interpolation
            % - replace land with the nearest point data -
            % find land/water cells
            idx_land=find(obc_mask < 0.5);
            idx_water=find(obc_mask > 0.5);
            % find water data at obc
            obc_lon1=obc_lon(idx_water);
            obc_lat1=obc_lat(idx_water);
            temi1=temi(idx_water);
            % find land data at obc
            obc_lon2=obc_lon(idx_land);
            obc_lat2=obc_lat(idx_land);
            % interpolate land with nearest data
            FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
            temi2=FN(obc_lon2,obc_lat2);
            % fill back to obc data
            temi(idx_land)=temi2; 
            zeta_east(tk,:)=temi;
        end
        % - western boundary -
        if obc_west
            obc_lon=lon_rho(:,1);
            obc_lat=lat_rho(:,1);
            obc_mask=mask_rho(:,1);
            temi=F(obc_lon,obc_lat); % interpolation
            % - replace land with the nearest point data -
            % find land/water cells
            idx_land=find(obc_mask < 0.5);
            idx_water=find(obc_mask > 0.5);
            % find water data at obc
            obc_lon1=obc_lon(idx_water);
            obc_lat1=obc_lat(idx_water);
            temi1=temi(idx_water);
            % find land data at obc
            obc_lon2=obc_lon(idx_land);
            obc_lat2=obc_lat(idx_land);
            % interpolate land with nearest data
            FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
            temi2=FN(obc_lon2,obc_lat2);
            % fill back to obc data
            temi(idx_land)=temi2; 
            zeta_west(tk,:)=temi;
        end
        % - northern boundary -
        if obc_north
            obc_lon=lon_rho(end,:)';
            obc_lat=lat_rho(end,:)';
            obc_mask=mask_rho(end,:)';
            temi=F(obc_lon,obc_lat); % interpolation
            % - replace land with the nearest point data -
            % find land/water cells
            idx_land=find(obc_mask < 0.5);
            idx_water=find(obc_mask > 0.5);
            % find water data at obc
            obc_lon1=obc_lon(idx_water);
            obc_lat1=obc_lat(idx_water);
            temi1=temi(idx_water);
            % find land data at obc
            obc_lon2=obc_lon(idx_land);
            obc_lat2=obc_lat(idx_land);
            % interpolate land with nearest data
            FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
            temi2=FN(obc_lon2,obc_lat2);
            % fill back to obc data
            temi(idx_land)=temi2; 
            zeta_north(tk,:)=temi;
        end
        % - southern boundary -
        if obc_south
            obc_lon=lon_rho(1,:)';
            obc_lat=lat_rho(1,:)';
            obc_mask=mask_rho(1,:)';
            temi=F(obc_lon,obc_lat); % interpolation
            % - replace land with the nearest point data -
            % find land/water cells
            idx_land=find(obc_mask < 0.5);
            idx_water=find(obc_mask > 0.5);
            % find water data at obc
            obc_lon1=obc_lon(idx_water);
            obc_lat1=obc_lat(idx_water);
            temi1=temi(idx_water);
            % find land data at obc
            obc_lon2=obc_lon(idx_land);
            obc_lat2=obc_lat(idx_land);
            % interpolate land with nearest data
            FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
            temi2=FN(obc_lon2,obc_lat2);
            % fill back to obc data
            temi(idx_land)=temi2; 
            zeta_south(tk,:)=temi;
        end     
%        
% -------- temperature and salinity interpolation ---------
        disp('           * Temperature and Salinity ...');
        % ----- horizontal interpolation -----
            
        for zk=1:nz
            % ===== salinity =====
            tem0=squeeze(salts(tk,zk,:,:));
            tem0=reshape(tem0,nx*ny,1);
            tem=[lon_data,lat_data,tem0]; % combine location and data together
            % remove NaN
            kk=find(isnan(sum(tem')));
            tem(kk,:)=[]; 
            tem=double(tem);            
            % - interpolaation -
            F = scatteredInterpolant(tem(:,1),tem(:,2),tem(:,3),'linear','nearest');
            % --- eastern boundary ---
            if obc_east
                obc_lon=lon_rho(:,end);
                obc_lat=lat_rho(:,end);
                obc_mask=mask_rho(:,end);
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obcer
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                salt_east0(zk,:)=temi;
            end
            % - western boundary -
            if obc_west
                obc_lon=lon_rho(:,1);
                obc_lat=lat_rho(:,1);
                obc_mask=mask_rho(:,1);
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                salt_west0(zk,:)=temi;
            end
            % - northern boundary -
            if obc_north
                obc_lon=lon_rho(end,:)';
                obc_lat=lat_rho(end,:)';
                obc_mask=mask_rho(end,:)';
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                salt_north0(zk,:)=temi;
            end
            % - southern boundary -
            if obc_south
                obc_lon=lon_rho(1,:)';
                obc_lat=lat_rho(1,:)';
                obc_mask=mask_rho(1,:)';
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                salt_south0(zk,:)=temi;
            end
        
            % ===== temperature =====
            tem0=squeeze(temps(tk,zk,:,:));
            tem0=reshape(tem0,nx*ny,1);
            tem=[lon_data,lat_data,tem0]; % combine location and data together
            % remove NaN
            kk=find(isnan(sum(tem')));
            tem(kk,:)=[]; 
            tem=double(tem);
            % - interpolaation -
            F = scatteredInterpolant(tem(:,1),tem(:,2),tem(:,3),'linear','nearest');
            % - eastern boundary -
            if obc_east
                obc_lon=lon_rho(:,end);
                obc_lat=lat_rho(:,end);
                obc_mask=mask_rho(:,end);
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                temp_east0(zk,:)=temi;
            end
            % - western boundary -
            if obc_west
                obc_lon=lon_rho(:,1);
                obc_lat=lat_rho(:,1);
                obc_mask=mask_rho(:,1);
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                temp_west0(zk,:)=temi;
            end
            % - northern boundary -
            if obc_north
                obc_lon=lon_rho(end,:)';
                obc_lat=lat_rho(end,:)';
                obc_mask=mask_rho(end,:)';
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                temp_north0(zk,:)=temi;
            end
            % - southern boundary -
            if obc_south
                obc_lon=lon_rho(1,:)';
                obc_lat=lat_rho(1,:)';
                obc_mask=mask_rho(1,:)';
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                temp_south0(zk,:)=temi;
            end
        end
        %
        % -------- vertical interpolation --------
        % --- Eastern boundary ---
        if obc_east
            z0=squeeze(zr(:,:,end)); % depth at the boundary
            salt_east(tk,:,:)=zeros(N,Mp);
            temp_east(tk,:,:)=zeros(N,Mp);
            obc_mask=mask_rho(:,end);
            idx_water=find(obc_mask==1);
            nwater=length(idx_water);
            for mk=1:nwater
                ik=idx_water(mk);
                zi=z0(:,ik);        % model depth
                % === salinity ===
                tem1=salt_east0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                salt_east(tk,:,ik)=temi;
                % === temperature ===
                tem1=temp_east0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                temp_east(tk,:,ik)=temi;
            end
        end
    %
        % ------ Western boundary ------
        if obc_west
            z0=squeeze(zr(:,:,1)); % depth at the boundary
            salt_west(tk,:,:)=zeros(N,Mp);
            temp_west(tk,:,:)=zeros(N,Mp);
            obc_mask=mask_rho(:,1);
            idx_water=find(obc_mask==1);
            nwater=length(idx_water);
            for mk=1:nwater
                ik=idx_water(mk);
                zi=z0(:,ik);        % model depth
                % === salinity ===
                tem1=salt_west0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                salt_west(tk,:,ik)=temi;
                % === temperature ===
                tem1=temp_west0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                temp_west(tk,:,ik)=temi;
            end
        end
    
        % --- Northern boundary ---
        if obc_north
            z0=squeeze(zr(:,end,:)); % depth at the boundary
            salt_north(tk,:,:)=zeros(N,Lp);
            temp_north(tk,:,:)=zeros(N,Lp);
            obc_mask=mask_rho(end,:);
            idx_water=find(obc_mask==1);
            nwater=length(idx_water);
            for mk=1:nwater
                ik=idx_water(mk);
                zi=z0(:,ik);        % model depth
                % === salinity ===
                tem1=salt_north0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                salt_north(tk,:,ik)=temi;
                % === temperature ===
                tem1=temp_north0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                temp_north(tk,:,ik)=temi;
            end
        end
        
        % --- Southern boundary ---
        if obc_south
            z0=squeeze(zr(:,1,:)); % depth at the boundary
            salt_south(tk,:,:)=zeros(N,Lp);
            temp_south(tk,:,:)=zeros(N,Lp);
            obc_mask=mask_rho(1,:);
            idx_water=find(obc_mask==1);
            nwater=length(idx_water);
            for mk=1:nwater
                ik=idx_water(mk);
                zi=z0(:,ik);        % model depth
                % === salinity ===
                tem1=salt_south0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                salt_south(tk,:,ik)=temi;
                % === temperature ===
                tem1=temp_south0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                temp_south(tk,:,ik)=temi;
            end
        end
%        
% ---------------------- u  ---------------------
    %
        disp(['           * u ...']);
    % ------ horizontal interpolation ------
        for zk=1:nz
            tem0=squeeze(us(tk,zk,:,:));
            tem0=reshape(tem0,nx*ny,1);
            tem=[lon_data,lat_data,tem0]; % combine location and data together
            % remove NaN
            kk=find(isnan(sum(tem')));
            tem(kk,:)=[]; 
            tem=double(tem);
            % - interpolaation -
            F = scatteredInterpolant(tem(:,1),tem(:,2),tem(:,3),'linear','nearest');
            % - eastern boundary -
            if obc_east
                obc_lon=lon_u(:,end);
                obc_lat=lat_u(:,end);
                obc_mask=mask_u(:,end);
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                u_east0(zk,:)=temi;
            end
            % - western boundary -
            if obc_west
                obc_lon=lon_u(:,1);
                obc_lat=lat_u(:,1);
                obc_mask=mask_u(:,1);
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                u_west0(zk,:)=temi;
            end
            % - northern boundary -
            if obc_north
                obc_lon=lon_u(end,:)';
                obc_lat=lat_u(end,:)';
                obc_mask=mask_u(end,:)';
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                u_north0(zk,:)=temi;
            end
            % - southern boundary -
            if obc_south
                obc_lon=lon_u(1,:)';
                obc_lat=lat_u(1,:)';
                obc_mask=mask_u(1,:)';
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                u_south0(zk,:)=temi;
            end
        end    
%        
    % ---------- vertical interpolation ------------
        % --- eastern boundary ---
        if obc_east
            z0=squeeze(zu(:,:,end)); % depth at the boundary
            u_east(tk,:,:)=zeros(N,Mp);
            obc_mask=mask_u(:,end);
            idx_water=find(obc_mask==1);
            nwater=length(idx_water);
            for mk=1:nwater
                ik=idx_water(mk);
                zi=z0(:,ik);        % model depth
                tem1=u_east0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                u_east(tk,:,ik)=temi;
            end
        end
    %
        % --- western boundary ---
        if obc_west
            z0=squeeze(zu(:,:,1)); % depth at the boundary
            u_west(tk,:,:)=zeros(N,Mp);
            obc_mask=mask_u(:,1);
            idx_water=find(obc_mask==1);
            nwater=length(idx_water);
            for mk=1:nwater
                ik=idx_water(mk);
                zi=z0(:,ik);        % model depth
                tem1=u_west0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                u_west(tk,:,ik)=temi;
            end
        end
    %
        %
        % --- northern boundary ---
        if obc_north
            z0=squeeze(zu(:,end,:)); % depth at theboundary
            u_north(tk,:,:)=zeros(N,L);
            obc_mask=mask_u(end,:);
            idx_water=find(obc_mask==1);
            nwater=length(idx_water);
            for mk=1:nwater
                ik=idx_water(mk);
                zi=z0(:,ik);        % model depth
                tem1=u_north0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                u_north(tk,:,ik)=temi;
            end
        end
       %
       % --- southern boundary ---
       if obc_south
            z0=squeeze(zu(:,1,:)); % depth at the boundary
            u_south(tk,:,:)=zeros(N,L);
            obc_mask=mask_u(1,:);
            idx_water=find(obc_mask==1);
            nwater=length(idx_water);
            for mk=1:nwater
                ik=idx_water(mk);
                zi=z0(:,ik);        % model depth
                tem1=u_south0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                u_south(tk,:,ik)=temi;
            end
        end

% ---------------------- v  ---------------------
    %
        disp(['           * v ...']);
    % ------ horizontal interpolation ------
        for zk=1:nz
            tem0=squeeze(vs(tk,zk,:,:));
            tem0=reshape(tem0,nx*ny,1);
            tem=[lon_data,lat_data,tem0]; % combine location and data together
            % remove NaN
            kk=find(isnan(sum(tem')));
            tem(kk,:)=[]; 
            tem=double(tem);
            % - interpolaation -
            F = scatteredInterpolant(tem(:,1),tem(:,2),tem(:,3),'linear','nearest');
            % - eastern boundary -
            if obc_east
                obc_lon=lon_v(:,end);
                obc_lat=lat_v(:,end);
                obc_mask=mask_v(:,end);
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                v_east0(zk,:)=temi;
            end
            % - western boundary -
            if obc_west
                obc_lon=lon_v(:,1);
                obc_lat=lat_v(:,1);
                obc_mask=mask_v(:,1);
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                v_west0(zk,:)=temi;
            end
            % - northern boundary -
            if obc_north
                obc_lon=lon_v(end,:)';
                obc_lat=lat_v(end,:)';
                obc_mask=mask_v(end,:)';
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                v_north0(zk,:)=temi;
            end
            % - southern boundary -
            if obc_south
                obc_lon=lon_v(1,:)';
                obc_lat=lat_v(1,:)';
                obc_mask=mask_v(1,:)';
                temi=F(obc_lon,obc_lat); % interpolation
                % - replace land with the nearest point data -
                % find land/water cells
                idx_land=find(obc_mask < 0.5);
                idx_water=find(obc_mask > 0.5);
                % find water data at obc
                obc_lon1=obc_lon(idx_water);
                obc_lat1=obc_lat(idx_water);
                temi1=temi(idx_water);
                % find land data at obc
                obc_lon2=obc_lon(idx_land);
                obc_lat2=obc_lat(idx_land);
                % interpolate land with nearest data
                FN = scatteredInterpolant(obc_lon1,obc_lat1,temi1,'nearest','nearest');
                temi2=FN(obc_lon2,obc_lat2);
                % fill back to obc data
                temi(idx_land)=temi2; 
                v_south0(zk,:)=temi;
            end
        end    
%        
    % ---------- vertical interpolation ------------
        % --- eastern boundary ---
        if obc_east
            z0=squeeze(zv(:,:,end)); % depth at the boundary
            v_east(tk,:,:)=zeros(N,M);
           obc_mask=mask_v(:,end);
            idx_water=find(obc_mask==1);
            nwater=length(idx_water);
            for mk=1:nwater
                ik=idx_water(mk);
                zi=z0(:,ik);        % model depth
                tem1=v_east0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                v_east(tk,:,ik)=temi;
            end
        end
    %
        % --- western boundary ---
        if obc_west
            z0=squeeze(zv(:,:,1)); % depth at the boundary
            v_west(tk,:,:)=zeros(N,M);
            obc_mask=mask_v(:,1);
            idx_water=find(obc_mask==1);
            nwater=length(idx_water);
            for mk=1:nwater
                ik=idx_water(mk);
                zi=z0(:,ik);        % model depth
                tem1=v_west0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                v_west(tk,:,ik)=temi;
            end
        end
    %
        %
        % --- northern boundary ---
        if obc_north
            z0=squeeze(zv(:,end,:)); % depth at theboundary
            v_north(tk,:,:)=zeros(N,Lp);
            obc_mask=mask_v(end,:);
            idx_water=find(obc_mask==1);
            nwater=length(idx_water);
            for mk=1:nwater
                ik=idx_water(mk);
                zi=z0(:,ik);        % model depth
                tem1=v_north0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                v_north(tk,:,ik)=temi;
            end
        end
       %
       % --- southern boundary ---
       if obc_south
            z0=squeeze(zv(:,1,:)); % depth at the boundary
            v_south(tk,:,:)=zeros(N,Lp);
            obc_mask=mask_v(1,:);
            idx_water=find(obc_mask==1);
            nwater=length(idx_water);
            for mk=1:nwater
                ik=idx_water(mk);
                zi=z0(:,ik);        % model depth
                tem1=v_south0(:,ik);  % data
                tem1=flipud(tem1);  % change data up down
                if sum(isnan(tem1)) > nz-1 % for grids without data, use 0 
                    temi=0*ones(N,1);
                else
                    kk=find(isnan(tem1));
                    tem1(kk)=[];
                    zz=sd; zz(kk)=[];
                    temi=interp1(zz,tem1,zi,'linear'); % linear interpolation
                    kk0=find(isnan(temi));  % find NaN for lower layer
                    if length(kk0) > 0
                       % replace NaN with the nearst data 
                       dkk0=diff(kk0);
                       if max(dkk0) > 1  % NaN at both surface and bottom
                           kkk=find(dkk0==max(dkk0));
                           kk1=kk0(1:kkk); % surface layer 
                           kk2=kk0(kkk+1:end); % bottom layer
                           temi(kk1)=temi(max(kk1)+1)*ones(size(kk1)); % for near bottom
                           temi(kk2)=temi(min(kk2)-1)*ones(size(kk2)); % for near surface 
                       else
                           if min(kk0) == 1
                              temi(kk0)=temi(max(kk0)+1)*ones(size(kk0)); % for near bottom
                           else
                              temi(kk0)=temi(min(kk0)-1)*ones(size(kk0)); % for near surface 
                           end
                       end
                    end
                end   
                v_south(tk,:,ik)=temi;
            end
        end
    end

    
% -------------------------------------------------------------------------
% --------------- Rotation of velocity (from E-N to XI-ETA) ---------------
% -------------------------------------------------------------------------
        disp('      -> Rotating velocity ...');
%
    % ------ eastern boundary ------
    if obc_east
        u_eastr=0*u_east; % define rotated velocities
        v_eastr=0*v_east;
        % - angle -
        ang0=ang(:,end);
        ang1=repmat(ang0',[N,1]);
        nrecord=size(u_east,1); % total number of data record
        for tk=1:nrecord
            ang2(tk,:,:)=ang1;
        end
        % - rotation interior -
        uin=u_east(:,:,2:end-1);
        vin=0.5*(v_east(:,:,1:end-1)+v_east(:,:,2:end));
        ur=uin.*cos(ang2(:,:,2:end-1))+vin.*sin(ang2(:,:,2:end-1));
        vr=vin.*cos(ang2(:,:,2:end-1))-uin.*sin(ang2(:,:,2:end-1));
        vrr=0.5*(vr(:,:,1:end-1)+vr(:,:,2:end));
        u_eastr(:,:,2:end-1)=ur;
        v_eastr(:,:,2:end-1)=vrr;
        % - add boundary points -
        ur1=uin(:,:,1).*cos(ang2(:,:,1))+vin(:,:,1).*sin(ang2(:,:,1));
        vr1=vin(:,:,1).*cos(ang2(:,:,1))-uin(:,:,1).*sin(ang2(:,:,1));
        u_eastr(:,:,1)=ur1;
        v_eastr(:,:,1)=vr1;
        ur2=uin(:,:,end).*cos(ang2(:,:,end))+vin(:,:,end).*sin(ang2(:,:,end));
        vr2=vin(:,:,end).*cos(ang2(:,:,end))-uin(:,:,end).*sin(ang2(:,:,end));
        u_eastr(:,:,end)=ur2;
        v_eastr(:,:,end)=vr2;
    end
    clear ang0 ang1 ang2 ur vr vrr ur1 vr1 ur2 vr2 uin vin;

    % --- western boundary ---
    if obc_west
        u_westr=0*u_west;
        v_westr=0*v_west;
        % - angle -
        ang0=ang(:,1);
        ang1=repmat(ang0',[N,1]);
        nrecord=size(u_west,1); % total number of data record
        for tk=1:nrecord
            ang2(tk,:,:)=ang1;
        end
        % - rotation interior -
        uin=u_west(:,:,2:end-1);
        vin=0.5*(v_west(:,:,1:end-1)+v_west(:,:,2:end));
        ur=uin.*cos(ang2(:,:,2:end-1))+vin.*sin(ang2(:,:,2:end-1));
        vr=vin.*cos(ang2(:,:,2:end-1))-uin.*sin(ang2(:,:,2:end-1));
        vrr=0.5*(vr(:,:,1:end-1)+vr(:,:,2:end));
        u_westr(:,:,2:end-1)=ur;
        v_westr(:,:,2:end-1)=vrr;
        % - add boundary points -
        ur1=uin(:,:,1).*cos(ang2(:,:,1))+vin(:,:,1).*sin(ang2(:,:,1));
        vr1=vin(:,:,1).*cos(ang2(:,:,1))-uin(:,:,1).*sin(ang2(:,:,1));
        u_westr(:,:,1)=ur1;
        v_westr(:,:,1)=vr1;
        ur2=uin(:,:,end).*cos(ang2(:,:,end))+vin(:,:,end).*sin(ang2(:,:,end));
        vr2=vin(:,:,end).*cos(ang2(:,:,end))-uin(:,:,end).*sin(ang2(:,:,end));
        u_westr(:,:,end)=ur2;
        v_westr(:,:,end)=vr2;
    end
    clear ang0 ang1 ang2 ur vr vrr ur1 vr1 ur2 vr2 uin vin;

    % --- northern boundary ---
    if obc_north
        u_northr=0*u_north;
        v_northr=0*v_north;
        % - angle -
        ang0=ang(end,:);
        ang1=repmat(ang0,[N,1]);
        nrecord=size(u_north,1); % total number of data record
        for tk=1:nrecord
            ang2(tk,:,:)=ang1;
        end
        % - rotation interior -
        vin=v_north(:,:,2:end-1);
        uin=0.5*(u_north(:,:,1:end-1)+u_north(:,:,2:end));
        ur=uin.*cos(ang2(:,:,2:end-1))+vin.*sin(ang2(:,:,2:end-1));
        vr=vin.*cos(ang2(:,:,2:end-1))-uin.*sin(ang2(:,:,2:end-1));
        urr=0.5*(ur(:,:,1:end-1)+ur(:,:,2:end));
        u_northr(:,:,2:end-1)=urr;
        v_northr(:,:,2:end-1)=vr;
        % - add boundary points -
        ur1=uin(:,:,1).*cos(ang2(:,:,1))+vin(:,:,1).*sin(ang2(:,:,1));
        vr1=vin(:,:,1).*cos(ang2(:,:,1))-uin(:,:,1).*sin(ang2(:,:,1));
        u_northr(:,:,1)=ur1;
        v_northr(:,:,1)=vr1;
        ur2=uin(:,:,end).*cos(ang2(:,:,end))+vin(:,:,end).*sin(ang2(:,:,end));
        vr2=vin(:,:,end).*cos(ang2(:,:,end))-uin(:,:,end).*sin(ang2(:,:,end));
        u_northr(:,:,end)=ur2;
        v_northr(:,:,end)=vr2;
    end
    clear ang0 ang1 ang2 ur vr vrr ur1 vr1 ur2 vr2 uin vin;

    % --- southern boundary ---
    if obc_south
        u_southr=0*u_south;
        v_southr=0*v_south;
        % - angle -
        ang0=ang(1,:);
        ang1=repmat(ang0,[N,1]);
        nrecord=size(u_south,1); % total number of data record
        for tk=1:nrecord
            ang2(tk,:,:)=ang1;
        end
        % - rotation interior -
        vin=v_south(:,:,2:end-1);
        uin=0.5*(u_south(:,:,1:end-1)+u_south(:,:,2:end));
        ur=uin.*cos(ang2(:,:,2:end-1))+vin.*sin(ang2(:,:,2:end-1));
        vr=vin.*cos(ang2(:,:,2:end-1))-uin.*sin(ang2(:,:,2:end-1));
        urr=0.5*(ur(:,:,1:end-1)+ur(:,:,2:end));
        u_southr(:,:,2:end-1)=urr;
        v_southr(:,:,2:end-1)=vr;
        % - add boundary points -
        ur1=uin(:,:,1).*cos(ang2(:,:,1))+vin(:,:,1).*sin(ang2(:,:,1));
        vr1=vin(:,:,1).*cos(ang2(:,:,1))-uin(:,:,1).*sin(ang2(:,:,1));
        u_southr(:,:,1)=ur1;
        v_southr(:,:,1)=vr1;
        ur2=uin(:,:,end).*cos(ang2(:,:,end))+vin(:,:,end).*sin(ang2(:,:,end));
        vr2=vin(:,:,end).*cos(ang2(:,:,end))-uin(:,:,end).*sin(ang2(:,:,end));
        u_southr(:,:,end)=ur2;
        v_southr(:,:,end)=vr2;
    end
    clear ang0 ang1 ang2 ur vr vrr ur1 vr1 ur2 vr2 uin vin;

%
% -------------------------------------------------------------------------
% ------------------------------ Save data --------------------------------
% -------------------------------------------------------------------------
%
    % ------ depth-mean velocities ------
    % ---- eastern boundary ---
    if obc_east
        ua_eastr=[];
        ua=u_eastr;
        zz0=squeeze(zw(:,:,end));
        dz=zz0(2:end,:)-zz0(1:end-1,:);
        nrecord=size(u_east,1); % total number of data record
        for tk=1:nrecord
            tem=squeeze(ua(tk,:,:));
            tem=sum(tem.*dz);
            ua_eastr(tk,:)=tem./(zz0(end,:)-zz0(1,:));
        end
%
        va_eastr=[];
        ua=v_eastr;
        zz0=squeeze(zw(:,:,end));
        zz0=0.5*(zz0(:,1:end-1)+zz0(:,2:end));
        dz=zz0(2:end,:)-zz0(1:end-1,:);
        for tk=1:nrecord
            tem=squeeze(ua(tk,:,:));
            tem=sum(tem.*dz);
            va_eastr(tk,:)=tem./(zz0(end,:)-zz0(1,:));
        end
    end

    % --- western boundary ---
    if obc_west
        ua_westr=[];
        ua=u_westr;
        zz0=squeeze(zw(:,:,1));
        dz=zz0(2:end,:)-zz0(1:end-1,:);
        nrecord=size(u_west,1); % total number of data record
        for tk=1:nrecord
            tem=squeeze(ua(tk,:,:));
            tem=sum(tem.*dz);
            ua_westr(tk,:)=tem./(zz0(end,:)-zz0(1,:));
        end
%
        va_westr=[];
        ua=v_westr;
        zz0=squeeze(zw(:,:,1));
        zz0=0.5*(zz0(:,1:end-1)+zz0(:,2:end));
        dz=zz0(2:end,:)-zz0(1:end-1,:);
        for tk=1:nrecord
            tem=squeeze(ua(tk,:,:));
            tem=sum(tem.*dz);
            va_westr(tk,:)=tem./(zz0(end,:)-zz0(1,:));
        end
    end

    % ---- northern boundary ---
    if obc_north
        ua_northr=[];
        ua=u_northr;
        zz0=squeeze(zw(:,end,:));
        zz0=0.5*(zz0(:,1:end-1)+zz0(:,2:end));
        dz=zz0(2:end,:)-zz0(1:end-1,:);
        nrecord=size(u_north,1); % total number of data record
        for tk=1:nrecord
            tem=squeeze(ua(tk,:,:));
            tem=sum(tem.*dz);
            ua_northr(tk,:)=tem./(zz0(end,:)-zz0(1,:));
        end
%
        va_northr=[];
        ua=v_northr;
        zz0=squeeze(zw(:,end,:));
        dz=zz0(2:end,:)-zz0(1:end-1,:);
        for tk=1:nrecord
            tem=squeeze(ua(tk,:,:));
            tem=sum(tem.*dz);
            va_northr(tk,:)=tem./(zz0(end,:)-zz0(1,:));
        end
    end

    % ---- southern boundary ---
    if obc_south
        ua_southr=[];
        ua=u_southr;
        zz0=squeeze(zw(:,end,:));
        zz0=0.5*(zz0(:,1:end-1)+zz0(:,2:end));
        dz=zz0(2:end,:)-zz0(1:end-1,:);
        nrecord=size(u_south,1); % total number of data record
        for tk=1:nrecord
            tem=squeeze(ua(tk,:,:));
            tem=sum(tem.*dz);
            ua_southr(tk,:)=tem./(zz0(end,:)-zz0(1,:));
        end
%
        va_southr=[];
        ua=v_southr;
        zz0=squeeze(zw(:,end,:));
        dz=zz0(2:end,:)-zz0(1:end-1,:);
        for tk=1:nrecord
            tem=squeeze(ua(tk,:,:));
            tem=sum(tem.*dz);
            va_southr(tk,:)=tem./(zz0(end,:)-zz0(1,:));
        end
    end
%
    % --- time ---
    t_zeta=tt;
    t_temp=tt;
    t_salt=tt;
    t_v2d=tt;
    t_v3d=tt;
    % --- rename variables ---
    clear u_east u_west u_north u_south v_east v_west v_north v_south;
%
    if obc_east
       u_east=u_eastr;
       v_east=v_eastr;
       ua_east=ua_eastr;
       va_east=va_eastr;
    end
    
    if obc_west
       u_west=u_westr;
       v_west=v_westr;
       ua_west=ua_westr;
       va_west=va_westr;
    end
    
    if obc_north
       u_north=u_northr;
       v_north=v_northr;
       ua_north=ua_northr;
       va_north=va_northr;
    end
    
    if obc_south
       u_south=u_southr;
       v_south=v_southr;
       ua_south=ua_southr;
       va_south=va_southr;
    end
%
% ---- save data -----
%
    disp('      -> Saving data ...');
   
    outfile0=outfile(fk,:);
    save(outfile0,  'ts','t_zeta','t_temp','t_salt','t_v2d','t_v3d'); 

    zeta_west=permute(zeta_west,[2,1]);
    zeta_east=permute(zeta_east,[2,1]);
    zeta_north=permute(zeta_north,[2,1]);
    salt_west=permute(salt_west,[3,2,1]);
    salt_east=permute(salt_east,[3,2,1]);
    salt_north=permute(salt_north,[3,2,1]);
    temp_west=permute(temp_west,[3,2,1]);
    temp_east=permute(temp_east,[3,2,1]);
    temp_north=permute(temp_north,[3,2,1]);
    u_west=permute(u_west,[3,2,1]);
    u_east=permute(u_east,[3,2,1]);
    u_north=permute(u_north,[3,2,1]);
    v_west=permute(v_west,[3,2,1]);
    v_east=permute(v_east,[3,2,1]);
    v_north=permute(v_north,[3,2,1]);
    ua_west=permute(ua_west,[2,1]);
    ua_east=permute(ua_east,[2,1]);
    ua_north=permute(ua_north,[2,1]);
    va_west=permute(va_west,[2,1]);
    va_east=permute(va_east,[2,1]);
    va_north=permute(va_north,[2,1]);

    if obc_east
         save(outfile0,'zeta_east','salt_east','temp_east','u_east','v_east','ua_east','va_east','-append');
    end
    if obc_west
         save(outfile0,'zeta_west','salt_west','temp_west','u_west','v_west','ua_west','va_west','-append');
    end
    if obc_north
         save(outfile0,'zeta_north','salt_north','temp_north','u_north','v_north','ua_north','va_north','-append');
    end
    if obc_south
         save(outfile0,'zeta_south','salt_south','temp_south','u_south','v_south','ua_south','va_south','-append');
    end
end
disp('  Done !');





