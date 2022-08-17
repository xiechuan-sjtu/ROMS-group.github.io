%  Generate atmosphere forcing file for RAS model
%  Input:  atmosphere forcing variables in ERA5
%  Output: ROMS foring.nc file
%
%  Written by Xie Chuan
%  Date: 2021-11-02
%
clc; 
clear all; 
close all;
% ROMS model name
title='RAS';
% Set year and file part
Year = 1998; 
stage='p2';
% Time parameters
mybasedate = datenum(1998,1,1);   %set reference time in the model
if strcmp(stage,'p1')
    days_begin=datenum(Year,1,1);
    days_end=datenum(Year,4,30);
elseif strcmp(stage,'p2')
    days_begin=datenum(Year,5,1);
    days_end=datenum(Year,8,31);
elseif strcmp(stage,'p3')    
    days_begin=datenum(Year,9,1);
    days_end=datenum(Year,12,31);   
end
%  --- select focring varibles (by true/false)----
idUwind = true ;
idVwind = true ;
idTair = true ;
idQair = true ;
idPair = true ;
idcloud = true ;
idrain = true ;
idswrad = false ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------- Set the path-------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ROMS model grid file
ROMS_grdfiles_dir=['G:\Ross_amundsen_roms_model\grid_file\'];
grdfilename=[ROMS_grdfiles_dir,'RAS_grd_32layer_new.nc'];     
%ERA5 rawdata file path
ERA5_forcingvars_dir=['G:\Ross_amundsen_roms_model\atmosphere_forcing_file\A_ERA5_rawdata\'];
aicefilename=[ERA5_forcingvars_dir,'era5_aice_',num2str(Year),'.nc'];
cloudfilename=[ERA5_forcingvars_dir,'era5_cloud_',num2str(Year),'.nc'];
swradfilename=[ERA5_forcingvars_dir,'era5_swrad_',num2str(Year),'.nc'];
Pairfilename=[ERA5_forcingvars_dir,'era5_pair_',num2str(Year),'.nc'];
rainfilename=[ERA5_forcingvars_dir,'era5_rain_',num2str(Year),'.nc'];
Tempfilename=[ERA5_forcingvars_dir,'era5_temp_',num2str(Year),'.nc'];
windfilename=[ERA5_forcingvars_dir,'era5_wind_',num2str(Year),'.nc'];
matfiledir=['G:\Ross_amundsen_roms_model\atmosphere_forcing_file\B_RAS_model_forcing_interpolation_variables\'];
%------------------------------ ERA5 variables ----------------------------
if Year >=2003
    latitude=ncread(windfilename,'latitude');  
    longitude=ncread(windfilename,'longitude');
    [lon_era5,lat_era5]=meshgrid(longitude,latitude);
    lon_era5=double(lon_era5');
    lat_era5=double(lat_era5');
    time_era5=ncread(windfilename,'time');
else
    latitude=ncread(windfilename,'latitude');  
    longitude=ncread(windfilename,'longitude');
    longitude(401:1:1201,:)=[];
    [lon_era5,lat_era5]=meshgrid(longitude,latitude);
    lon_era5=double(lon_era5');
    lat_era5=double(lat_era5');
    time_era5=ncread(windfilename,'time');
end
% --------------------------- RAS Grid variables --------------------------
ang_ras=ncread(grdfilename,'angle');
lon_rho=ncread(grdfilename,'lon_rho');
lat_rho=ncread(grdfilename,'lat_rho');
%longitude convert to positive value
if Year >= 2003
    for i=1:1000
        for j=1:1000
             if lon_rho(i,j)<0
                 lon_rho(i,j)=lon_rho(i,j)+360;
             end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------- read atmospheric forcing field variables ---------------
% ------------------- interpolated into the RAS model grid ----------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------------------------time frequency------------------------------
%--------------------------------------------------------------------------
dt = 1/8;  %time frequency (unit:days)  here： 3 hours
interval = dt*24-1;
time_end = (days_end-days_begin+1)./dt;
reftime = days_begin-mybasedate;
%--------------------------------------------------------------------------
%---------------------------------Uwind------------------------------------
if idUwind
    Uwind=ncread(windfilename,'u10');
    if Year <= 2002
        Uwind(401:1:1201,:,:)=[];
    end
    for i=1:time_end
        mark_index=(i-1)*(24*dt)+1+(days_begin-datenum(Year,1,1))*24;
        time_now=(i-1)*dt+reftime;
        wind_time(i,1)=time_now;
        disp(['Interpolate Uwind in RAS model at time index  : ',num2str(time_now)]);
        disp(['Interpolate Uwind from ERA5 dataset at index  : ',num2str(mark_index)]);
        Uwind_tmp=mean(Uwind(:,:,mark_index:mark_index+interval),3);
        F=griddata(lon_era5,lat_era5,double(Uwind_tmp),lon_rho,lat_rho,'cubic');
        Uwind_RAS(:,:,i)=fillmissing(F,'nearest');  
        clear F; 
    end
    clear Uwind
    matfilename=[matfiledir,'RAS_Uwind_',num2str(Year),stage,'.mat'];
    save(matfilename,'Uwind_RAS','wind_time');
    clear Uwind_RAS
end
%--------------------------------------------------------------------------
%---------------------------------Vwind------------------------------------
if idVwind
    Vwind=ncread(windfilename,'v10');
    if Year <= 2002
       Vwind(401:1:1201,:,:)=[];
    end
    for i=1:time_end
        mark_index=(i-1)*(24*dt)+1+(days_begin-datenum(Year,1,1))*24;
        time_now=(i-1)*dt+reftime;
        wind_time(i,1)=time_now;
        disp(['Interpolate Vwind in RAS model at time index  : ',num2str(time_now)]);
        disp(['Interpolate Vwind from ERA5 dataset at index  : ',num2str(mark_index)]);
        Vwind_tmp=mean(Vwind(:,:,mark_index:mark_index+interval),3);
        F=griddata(lon_era5,lat_era5,double(Vwind_tmp),lon_rho,lat_rho,'cubic');
        Vwind_RAS(:,:,i)=fillmissing(F,'nearest');  
        clear F;
    end
    clear Vwind
    matfilename=[matfiledir,'RAS_Vwind_',num2str(Year),stage,'.mat'];
    save(matfilename,'Vwind_RAS','wind_time');
    clear Vwind_RAS
end
%--------------------------------------------------------------------------
%-----------------------------------Tair-----------------------------------
if idTair
    Tair=ncread(Tempfilename,'t2m')-273.15;   % Kelvin to Celsius
    if Year <= 2002
       Tair(401:1:1201,:,:)=[];
    end
    for i=1:time_end
        mark_index=(i-1)*(24*dt)+1+(days_begin-datenum(Year,1,1))*24;
        time_now=(i-1)*dt+reftime;
        tair_time(i,1)=time_now;
        disp(['Interpolate Tair in RAS model at time index  : ',num2str(time_now)]);
        disp(['Interpolate Tair from ERA5 dataset at index  : ',num2str(mark_index)]);
        Tair_tmp=mean(Tair(:,:,mark_index:mark_index+interval),3);
        F=griddata(lon_era5,lat_era5,double(Tair_tmp),lon_rho,lat_rho,'cubic');
        Tair_RAS(:,:,i)=fillmissing(F,'nearest');  
        clear F;
    end
    clear Tair
    matfilename=[matfiledir,'RAS_Tair_',num2str(Year),stage,'.mat'];
    save(matfilename,'Tair_RAS','tair_time');
    clear Tair_RAS
end
%-------------------------End time index----------------------------------
mark_index_end=(time_end-1)*(24*dt)+1+(days_begin-datenum(Year,1,1))*24;
filetime_end=(time_end-1)*dt+reftime;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
dt=2;  %time frequency (unit:days)  here: 2 day
interval=dt*24-1;
time_end=(days_end-days_begin+1)./dt;
reftime=days_begin-mybasedate;
%--------------------------------------------------------------------------
%----------------------------------Qair------------------------------------
if idQair
    tsur  = ncread(Tempfilename, 't2m')- 273.15;     % 2m temperature
    if Year <= 2002
       tsur(401:1:1201,:,:)=[];
    end
    tdew  = ncread(Tempfilename, 'd2m')- 273.15;     % 2m dewpoint
    if Year <= 2002
       tdew(401:1:1201,:,:)=[];
    end
    E     = 6.11 .* 10.0 .^ (7.5 .* tdew ./ (237.7 + tdew));
    clear tdew
    Es    = 6.11 .* 10.0 .^ (7.5 .* tsur ./ (237.7 + tsur));
    clear tsur
    Qair  = 100.0 .* (E ./ Es); 
    clear E Es
    for i=1:time_end
        mark_index=(i-1)*(24*dt)+1+(days_begin-datenum(Year,1,1))*24;
        time_now=(i-1)*dt+reftime;
        qair_time(i,1)=time_now;
        disp(['Interpolate Qair in RAS model at time index  : ',num2str(time_now)]);
        disp(['Interpolate Qair from ERA5 dataset at index  : ',num2str(mark_index)]);
        Qair_tmp=mean(Qair(:,:,mark_index:mark_index+interval),3);
        F=griddata(lon_era5,lat_era5,double(Qair_tmp),lon_rho,lat_rho,'cubic');
        Qair_RAS(:,:,i)=fillmissing(F,'nearest');  
        clear F;
    end
    [dim,~]=size(qair_time);
    disp(['Interpolate Qair in RAS model at time index  : ',num2str(filetime_end)]);
    disp(['Interpolate Qair from ERA5 dataset at index  : ',num2str(mark_index_end)]);
    qair_time(dim+1,1)=filetime_end;
    Qair_RAS(:,:,dim+1)=Qair_RAS(:,:,dim);
    % F=griddata(lon_era5,lat_era5,double(Qair(:,:,mark_index_end)),lon_rho,lat_rho,'nearest');
    % Qair_RAS(:,:,dim+1)=fillmissing(F,'nearest');
    clear F Qair
    matfilename=[matfiledir,'RAS_Qair_',num2str(Year),stage,'.mat'];
    save(matfilename,'Qair_RAS','qair_time');
    clear Qair_RAS
end
%--------------------------------------------------------------------------
%----------------------------------cloud-----------------------------------
if idcloud
    cloud=ncread(cloudfilename,'tcc');
    if Year <= 2002
       cloud(401:1:1201,:,:)=[];
    end
    for i=1:time_end
        mark_index=(i-1)*(24*dt)+1+(days_begin-datenum(Year,1,1))*24;
        time_now=(i-1)*dt+reftime;
        cloud_time(i,1)=time_now;
        disp(['Interpolate cloud in RAS model at time index  : ',num2str(time_now)]);
        disp(['Interpolate cloud from ERA5 dataset at index  : ',num2str(mark_index)]);
        cloud_tmp=mean(cloud(:,:,mark_index:mark_index+interval),3);
        F=griddata(lon_era5,lat_era5,double(cloud_tmp),lon_rho,lat_rho,'cubic');
        cloud_RAS(:,:,i)=fillmissing(F,'nearest');  
        clear F;
    end
    [dim,~]=size(cloud_time);
    disp(['Interpolate cloud in RAS model at time index  : ',num2str(filetime_end)]);
    disp(['Interpolate cloud from ERA5 dataset at index  : ',num2str(mark_index_end)]);
    cloud_time(dim+1,1)=filetime_end;
    cloud_RAS(:,:,dim+1)=cloud_RAS(:,:,dim);
    % F=griddata(lon_era5,lat_era5,double(cloud(:,:,mark_index_end)),lon_rho,lat_rho,'nearest');
    % cloud_RAS(:,:,dim+1)=fillmissing(F,'nearest');
    clear F cloud
    matfilename=[matfiledir,'RAS_cloud_',num2str(Year),stage,'.mat'];
    save(matfilename,'cloud_RAS','cloud_time');
    clear cloud_RAS
end
%--------------------------------------------------------------------------
%----------------------------------swrad-----------------------------------
if idswrad
    swrad=ncread(swradfilename,'msdwswrf');
    if Year <= 2002
       swrad(401:1:1201,:,:)=[];
    end
    for i=1:time_end
        mark_index=(i-1)*(24*dt)+1+(days_begin-datenum(Year,1,1))*24;
        time_now=(i-1)*dt+reftime;
        srf_time(i,1)=time_now;
        disp(['Interpolate swrad in RAS model at time index  : ',num2str(time_now)]);
        disp(['Interpolate swrad from ERA5 dataset at index  : ',num2str(mark_index)]);
        swrad_tmp=mean(swrad(:,:,mark_index:mark_index+interval),3);
        F=griddata(lon_era5,lat_era5,double(swrad_tmp),lon_rho,lat_rho,'cubic');
        swrad_RAS(:,:,i)=fillmissing(F,'nearest');  
        clear F;
    end
    [dim,~]=size(srf_time);
    disp(['Interpolate swrad in RAS model at time index  : ',num2str(filetime_end)]);
    disp(['Interpolate swrad from ERA5 dataset at index  : ',num2str(mark_index_end)]);
    srf_time(dim+1,1)=filetime_end;
    swrad_RAS(:,:,dim+1)=swrad_RAS(:,:,dim);
    % F=griddata(lon_era5,lat_era5,double(swrad(:,:,mark_index_end)),lon_rho,lat_rho,'nearest');
    % swrad_RAS(:,:,dim+1)=fillmissing(F,'nearest');
    clear F swrad
    matfilename=[matfiledir,'RAS_swrad_',num2str(Year),stage,'.mat'];
    save(matfilename,'swrad_RAS','srf_time');
    clear swrad_RAS
end
%--------------------------------------------------------------------------
%-----------------------------------Pair-----------------------------------
if idPair
    Pair=ncread(Pairfilename,'msl')./100;
    if Year <= 2002
       Pair(401:1:1201,:,:)=[];
    end
    for i=1:time_end
        mark_index=(i-1)*(24*dt)+1+(days_begin-datenum(Year,1,1))*24;
        time_now=(i-1)*dt+reftime;
        pair_time(i,1)=time_now;
        disp(['Interpolate Pair in RAS model at time index  : ',num2str(time_now)]);
        disp(['Interpolate Pair from ERA5 dataset at index  : ',num2str(mark_index)]);
        Pair_tmp=mean(Pair(:,:,mark_index:mark_index+interval),3);
        F=griddata(lon_era5,lat_era5,double(Pair_tmp),lon_rho,lat_rho,'cubic');
        Pair_RAS(:,:,i)=fillmissing(F,'nearest');  
        clear F;
    end
    [dim,~]=size(pair_time);
    disp(['Interpolate Pair in RAS model at time index  : ',num2str(filetime_end)]);
    disp(['Interpolate Pair from ERA5 dataset at index  : ',num2str(mark_index_end)]);
    pair_time(dim+1,1)=filetime_end;
    Pair_RAS(:,:,dim+1)=Pair_RAS(:,:,dim);
    % F=griddata(lon_era5,lat_era5,double(Pair(:,:,mark_index_end)),lon_rho,lat_rho,'nearest');
    % Pair_RAS(:,:,dim+1)=fillmissing(F,'nearest');
    clear F Pair
    matfilename=[matfiledir,'RAS_Pair_',num2str(Year),stage,'.mat'];
    save(matfilename,'Pair_RAS','pair_time');
    clear Pair_RAS
end
%--------------------------------------------------------------------------
%-----------------------------------rain-----------------------------------
if idrain
    rain=ncread(rainfilename,'mtpr');
    if Year <= 2002
       rain(401:1:1201,:,:)=[];
    end
    for i=1:time_end
        mark_index=(i-1)*(24*dt)+1+(days_begin-datenum(Year,1,1))*24;
        time_now=(i-1)*dt+reftime;
        rain_time(i,1)=time_now;
        disp(['Interpolate rain in RAS model at time index  : ',num2str(time_now)]);
        disp(['Interpolate rain from ERA5 dataset at index  : ',num2str(mark_index)]);
        rain_tmp=mean(rain(:,:,mark_index:mark_index+interval),3);
        F=griddata(lon_era5,lat_era5,double(rain_tmp),lon_rho,lat_rho,'cubic');
        rain_RAS(:,:,i)=fillmissing(F,'nearest');  
        clear F;
    end
    [dim,~]=size(rain_time);
    disp(['Interpolate rain in RAS model at time index  : ',num2str(filetime_end)]);
    disp(['Interpolate rain from ERA5 dataset at index  : ',num2str(mark_index_end)]);
    rain_time(dim+1,1)=filetime_end;
    rain_RAS(:,:,dim+1)=rain_RAS(:,:,dim);
    % F=griddata(lon_era5,lat_era5,double(rain(:,:,mark_index_end)),lon_rho,lat_rho,'nearest');
    % rain_RAS(:,:,dim+nearest1)=fillmissing(F,'');
    clear F rain
    matfilename=[matfiledir,'RAS_rain_',num2str(Year),stage,'.mat'];
    save(matfilename,'rain_RAS','rain_time');
    clear rain_RAS
end

%% wind ratation
ROMS_grdfiles_dir=['G:\Ross_amundsen_roms_model\grid_file\'];
grdfilename=[ROMS_grdfiles_dir,'RAS_grd_32layer_new.nc'];     
ang_ras=ncread(grdfilename,'angle');
stage=['p1';'p2';'p3'];
Year=1998;  
for i=1:3
    tic
    load(['RAS_Uwind_',num2str(Year),stage(i,:),'.mat']);
    load(['RAS_Vwind_',num2str(Year),stage(i,:),'.mat']);
    [time_end,~]=size(wind_time);
    %-----wind rotation (xi-eta)----
    for indext=1:time_end
        disp(indext)
        cff1=Uwind_RAS(:,:,indext);
        cff2=Vwind_RAS(:,:,indext);
        Uwindxi_RAS(:,:,indext)=cff1.*cos(ang_ras)+cff2.*sin(ang_ras);
        Vwindeta_RAS(:,:,indext)=cff2.*cos(ang_ras)-cff1.*sin(ang_ras);
    end
    clear Uwind_RAS Vwind_RAS
    save(['RAS_Uwindxi_',num2str(Year),stage(i,:),'.mat'],'Uwindxi_RAS','wind_time');
    save(['RAS_Vwindeta_',num2str(Year),stage(i,:),'.mat'],'Vwindeta_RAS','wind_time');
    clear Uwindxi_RAS Vwindeta_RAS
    disp('Done!')
    toc
end
%%
% Write netcdf files 格式 NETCDF4
% clc;clear all;
clc;
clear all;

Year=2014;
stage='p3';
%ROMS网格文件路径
ROMS_grdfiles_dir=['G:\Ross_amundsen_roms_model\grid_file\'];
grdfilename=[ROMS_grdfiles_dir,'RAS_grd_32layer_new.nc'];     
ROMS_forcingfiles_dir=['G:\Ross_amundsen_roms_model\atmosphere_forcing_file\E_High_time_resolution_forcing_files\'];
matfiledir=['G:\Ross_amundsen_roms_model\atmosphere_forcing_file\D_RAS_model_forcing_interpolation_variables\'];
ncfilename=[ROMS_forcingfiles_dir,'RAS_forcing_',num2str(Year),stage,'_era5.nc'];
[londim,~]=size(ncread(grdfilename,'lon_rho'));
[~,latdim]=size(ncread(grdfilename,'lat_rho'));
% ---------------------- creat file ----------------------
disp('  Creating NC file...')
% -----------------------load data -----------------------
load([matfiledir,'RAS_Uwind_',num2str(Year),stage,'.mat']);clear Uwind_RAS;
% load([matfiledir,'RAS_Tair_',num2str(Year),stage,'.mat']);clear Tair_RAS;
load([matfiledir,'RAS_Qair_',num2str(Year),stage,'.mat']);
load([matfiledir,'RAS_cloud_',num2str(Year),stage,'.mat']);
load([matfiledir,'RAS_swrad_',num2str(Year),stage,'.mat']);
load([matfiledir,'RAS_Pair_',num2str(Year),stage,'.mat']);
load([matfiledir,'RAS_rain_',num2str(Year),stage,'.mat']);
% ------------------- define dimensions -------------------
disp('      => Defining Dimensions...')
dim_xi_rho=londim;
dim_eta_rho=latdim;
[dim_wind_time,~]=size(wind_time);
[dim_tair_time,~]=size(wind_time);
[dim_pair_time,~]=size(pair_time);
[dim_qair_time,~]=size(qair_time);
[dim_cloud_time,~]=size(cloud_time);
[dim_srf_time,~]=size(srf_time);
[dim_rain_time,~]=size(rain_time);
% ------------------- define Variables and attributes -------------------
disp('      => Defining Variables and Attributes...')
% --------------------------------------------------------------------------
% ------------------------------- wind_time -------------------------------- 
nccreate(ncfilename,'wind_time','Dimensions',{'wind_time',dim_wind_time},'FillValue','disable');
ncwriteatt(ncfilename,'wind_time','long_name','surface wind time');
ncwriteatt(ncfilename,'wind_time','units','day');
ncwriteatt(ncfilename,'wind_time','field','wind_time, scalar, series');
% ------------------------------- tair_time -------------------------------- 
nccreate(ncfilename,'tair_time','Dimensions',{'tair_time',dim_tair_time},'FillValue','disable');
ncwriteatt(ncfilename,'tair_time','long_name','surface air temperature time');
ncwriteatt(ncfilename,'tair_time','units','day');
ncwriteatt(ncfilename,'tair_time','field','tair_time, scalar, series');
% ------------------------------- pair_time -------------------------------- 
nccreate(ncfilename,'pair_time','Dimensions',{'pair_time',dim_pair_time},'FillValue','disable');
ncwriteatt(ncfilename,'pair_time','long_name','surface air pressure time');
ncwriteatt(ncfilename,'pair_time','units','day');
ncwriteatt(ncfilename,'pair_time','field','pair_time, scalar, series');
% ------------------------------- qair_time -------------------------------- 
nccreate(ncfilename,'qair_time','Dimensions',{'qair_time',dim_qair_time},'FillValue','disable');
ncwriteatt(ncfilename,'qair_time','long_name','surface air relative humidity time');
ncwriteatt(ncfilename,'qair_time','units','day');
ncwriteatt(ncfilename,'qair_time','field','qair_time, scalar, series');
% ------------------------------ cloud_time ------------------------------- 
nccreate(ncfilename,'cloud_time','Dimensions',{'cloud_time',dim_cloud_time},'FillValue','disable');
ncwriteatt(ncfilename,'cloud_time','long_name','cloud fraction time');
ncwriteatt(ncfilename,'cloud_time','units','day');
ncwriteatt(ncfilename,'cloud_time','field','cloud_time, scalar, series');
% ------------------------------- srf_time -------------------------------- 
nccreate(ncfilename,'srf_time','Dimensions',{'srf_time',dim_srf_time},'FillValue','disable');
ncwriteatt(ncfilename,'srf_time','long_name','solar shortwave radiation flux time');
ncwriteatt(ncfilename,'srf_time','units','day');
ncwriteatt(ncfilename,'srf_time','field','srf_time, scalar, series');
% ------------------------------- rain_time -------------------------------- 
nccreate(ncfilename,'rain_time','Dimensions',{'rain_time',dim_rain_time},'FillValue','disable');
ncwriteatt(ncfilename,'rain_time','long_name','rain fall rate time');
ncwriteatt(ncfilename,'rain_time','units','day');
ncwriteatt(ncfilename,'rain_time','field','rain_time, scalar, series');
% ------------------------------- Uwind -------------------------------- 
nccreate(ncfilename,'Uwind','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'wind_time',dim_wind_time},'FillValue','disable');
ncwriteatt(ncfilename,'Uwind','long_name','surface u-wind component');
ncwriteatt(ncfilename,'Uwind','units','meter second-1');
ncwriteatt(ncfilename,'Uwind','field','u-wind, scalar, series');
ncwriteatt(ncfilename,'Uwind','time','wind_time');
% ------------------------------- Vwind -------------------------------- 
nccreate(ncfilename,'Vwind','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'wind_time',dim_wind_time},'FillValue','disable');
ncwriteatt(ncfilename,'Vwind','long_name','surface v-wind component');
ncwriteatt(ncfilename,'Vwind','units','meter second-1');
ncwriteatt(ncfilename,'Vwind','field','v-wind, scalar, series');
ncwriteatt(ncfilename,'Vwind','time','wind_time');
% ------------------------------- Tair -------------------------------- 
nccreate(ncfilename,'Tair','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'tair_time',dim_tair_time},'FillValue','disable');
ncwriteatt(ncfilename,'Tair','long_name','surface air temperature');
ncwriteatt(ncfilename,'Tair','units','Celsius');
ncwriteatt(ncfilename,'Tair','field','Tair, scalar, series');
ncwriteatt(ncfilename,'Tair','time','tair_time');
% ------------------------------- Pair -------------------------------- 
nccreate(ncfilename,'Pair','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'pair_time',dim_pair_time},'FillValue','disable');
ncwriteatt(ncfilename,'Pair','long_name','surface air pressure');
ncwriteatt(ncfilename,'Pair','units','milibar');
ncwriteatt(ncfilename,'Pair','field','Pair, scalar, series');
ncwriteatt(ncfilename,'Pair','time','pair_time');
% ------------------------------- Qair -------------------------------- 
nccreate(ncfilename,'Qair','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'qair_time',dim_qair_time},'FillValue','disable');
ncwriteatt(ncfilename,'Qair','long_name','surface air relative humidity');
ncwriteatt(ncfilename,'Qair','units','percentage');
ncwriteatt(ncfilename,'Qair','field','Qair, scalar, series');
ncwriteatt(ncfilename,'Qair','time','qair_time');
% ------------------------------- cloud -------------------------------- 
nccreate(ncfilename,'cloud','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'cloud_time',dim_cloud_time},'FillValue','disable');
ncwriteatt(ncfilename,'cloud','long_name','cloud fraction');
ncwriteatt(ncfilename,'cloud','units','nondimensional');
ncwriteatt(ncfilename,'cloud','field','cloud, scalar, series');
ncwriteatt(ncfilename,'cloud','time','cloud_time');
% ------------------------------- swrad -------------------------------- 
nccreate(ncfilename,'swrad','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'srf_time',dim_srf_time},'FillValue','disable');
ncwriteatt(ncfilename,'swrad','long_name','solar shortwave radiation flux');
ncwriteatt(ncfilename,'swrad','units','watt meter-2');
ncwriteatt(ncfilename,'swrad','field','shortwave radiation, scalar, series');
ncwriteatt(ncfilename,'swrad','time','srf_time');
% ------------------------------- rain -------------------------------- 
nccreate(ncfilename,'rain','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'rain_time',dim_rain_time},'FillValue','disable');
ncwriteatt(ncfilename,'rain','long_name','rain fall rate');
ncwriteatt(ncfilename,'rain','units','kilogram meter-2 second-1');
ncwriteatt(ncfilename,'rain','field','rain, scalar, series');
ncwriteatt(ncfilename,'rain','time','rain_time');
%--------------------------------------------------------------------------
disp('      => Loading Forcing data ...');
%--------------------------------------------------------------------------
% ---------------- write variables -----------------
disp('      => Filling Variables...')

load([matfiledir,'RAS_Uwind_',num2str(Year),stage,'.mat']);
ncwrite(ncfilename,'wind_time',wind_time);
ncwrite(ncfilename,'Uwind',Uwind_RAS);clear Uwind_RAS;
load([matfiledir,'RAS_Vwind_',num2str(Year),stage,'.mat']);
ncwrite(ncfilename,'Vwind',Vwind_RAS);clear Vwind_RAS;
load([matfiledir,'RAS_Tair_',num2str(Year),stage,'.mat']);
ncwrite(ncfilename,'tair_time',tair_time);
ncwrite(ncfilename,'Tair',Tair_RAS);clear Tair_RAS;
ncwrite(ncfilename,'qair_time',qair_time);
ncwrite(ncfilename,'Qair',Qair_RAS);clear Qair_RAS;
ncwrite(ncfilename,'cloud_time',cloud_time);
ncwrite(ncfilename,'cloud',cloud_RAS);clear cloud_RAS;
ncwrite(ncfilename,'srf_time',srf_time);
ncwrite(ncfilename,'swrad',swrad_RAS);clear swrad_RAS;
ncwrite(ncfilename,'pair_time',pair_time);
ncwrite(ncfilename,'Pair',Pair_RAS);clear Pair_RAS;
ncwrite(ncfilename,'rain_time',rain_time);
ncwrite(ncfilename,'rain',rain_RAS);clear rain_RAS;

%----------------- Write Global Attributes ----------------
disp('      => Creating Global Attributes...')
ncwriteatt(ncfilename,'/','History','FORCING file for RAS model');
ncwriteatt(ncfilename,'/','Create_time',datestr(now));
ncwriteatt(ncfilename,'/','Type','ERA5 FORCING file');
ncwriteatt(ncfilename,'/','Time','days (since 2003-01-01 00:00:00)');
ncwriteatt(ncfilename,'/','Out_file',ncfilename);
ncwriteatt(ncfilename,'/','Grd_file',grdfilename);

disp('  Done !');

