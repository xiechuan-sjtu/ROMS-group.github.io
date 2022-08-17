% Generate initial condition for ROMS
% The model starts from rest
% Usage
%  Input:  initial temperatures and salinities,
%  Output: ROMS ini.nc file
%
% Written by Xie Chuan,10-11-2021
% Relevant functions are stored in the utitilty fold under Matlab/toolbox
%
clear all; close all; clc;
% -------------------------------------------------------------------------
% ------------------------- user settings ---------------------------------
% -------------------------------------------------------------------------
%
title='RAS';
%
% Common parameters
%
% romstools_param
ROMS_grdfiles_dir=['F:\Ross_amundsen_roms_model\grid_file\'];
ROMS_inifiles_dir=['F:\Ross_amundsen_roms_model\initialization_file\'];
grdfilename=[ROMS_grdfiles_dir,'RAS_grd_32layer.nc'];
inifilename=[ROMS_inifiles_dir,'RAS_ini.nc'];
%
% --- S-coordinate setup ---
%
N =32;           % No. of vertical layers
Vtransform = 2;   % Vertical transformation equation: 1, original; 2, new
Vstretching = 4;  % Vstretching   Vertical stretching function:
%                    Vstretching = 1,  original (Song and Haidvogel, 1994)
%                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS, 2005)
%                    Vstretching = 3,  R. Geyer BBL refinement
%                    Vstretching = 4,  A. Shchepetkin (UCLA-ROMS, 2010)
theta_s = ncread(grdfilename,'theta_s');    % S-coordinate surface control parameter (scalar)
theta_b = ncread(grdfilename,'theta_b');    % S-coordinate bottom control parameter (scalar)
Tcline = 200; 
hc = 200;         % Width (m) of surface or bottom boundary layer in which
                  %  higher vertical resolution is required during
                  %  stretching (scalar)
%
%-----Read the grid file-------
%
spherical = ncread(grdfilename,'spherical');
s_w = ncread(grdfilename,'s_w');
s_rho = ncread(grdfilename,'s_rho');
Cs_r = ncread(grdfilename,'Cs_r');
Cs_w = ncread(grdfilename,'Cs_w');
h = ncread(grdfilename,'h');
lon_rho = ncread(grdfilename,'lon_rho');
lat_rho = ncread(grdfilename,'lat_rho');
lon_u = ncread(grdfilename,'lon_u');
lat_u = ncread(grdfilename,'lat_u');
lon_v = ncread(grdfilename,'lon_v');
lat_v = ncread(grdfilename,'lat_v');
[Mp,Lp]=size(h);
L=Lp-1;
M=Mp-1;
Np=N+1;
%
% *************************************************************************
% ********************** No Changes needed below **************************
% *************************************************************************
%
% ------------------------ Generate initial field -------------------------
% 
disp('  Generating initial data ...');
% initial time
% ocean_time_ini=datenum(yr_ini,mon_ini,day_ini)-datenum(yr_ref,mon_ref,day_ref);
ocean_time_ini=0;
%
% ---- interpolation ----
disp('  Generating initial fields ...');
% -------------------------------------------------------------------------
% ------------------------- Create netcdf file ----------------------------
% -------------------------------------------------------------------------
% --- creat file ---
disp('  Creating NC file...')
% --- define dimensions ---
disp('    => Defining Dimensions...')
% %
dim_xi_rho = Lp;
dim_xi_u = L;
dim_xi_v = Lp;
dim_eta_rho = Mp;
dim_eta_u = Mp;
dim_eta_v = M; 
dim_s_rho = N;
dim_s_w = Np;
dim_tracer = 5;
dim_ocean_time = 1;
dim_one = 1;
% 
% --- define Variables and attributes ---
disp('    => Defining Variables and Attributes...')
% %--------------------------------------------------------------------------
nccreate(inifilename,'spherical','Dimensions',{'one',dim_one},'FillValue','disable');
ncwriteatt(inifilename,'spherical','long_name','grid type logical switch')
ncwriteatt(inifilename,'spherical','flag_values',[0  1])
ncwriteatt(inifilename,'spherical','flag_meanings','Cartesian spherical')
% %--------------------------------------------------------------------------
nccreate(inifilename,'Vtransform','Dimensions',{'one',dim_one},'FillValue','disable');
ncwriteatt(inifilename,'Vtransform','long_name','vertical terrain-following transformation equation')
% %--------------------------------------------------------------------------
nccreate(inifilename,'Vstretching','Dimensions',{'one',dim_one},'FillValue','disable');
ncwriteatt(inifilename,'Vstretching','long_name','vertical terrain-following stretching function')
% %--------------------------------------------------------------------------
nccreate(inifilename,'theta_s','Dimensions',{'one',dim_one},'FillValue','disable');
ncwriteatt(inifilename,'theta_s','long_name','S-coordinate surface control parameter')
% %--------------------------------------------------------------------------
nccreate(inifilename,'theta_b','Dimensions',{'one',dim_one},'FillValue','disable');
ncwriteatt(inifilename,'theta_b','long_name','S-coordinate bottom control parameter')
% %--------------------------------------------------------------------------
nccreate(inifilename,'Tcline','Dimensions',{'one',dim_one},'FillValue','disable');
ncwriteatt(inifilename,'Tcline','long_name','S-coordinate surface/bottom layer width')
ncwriteatt(inifilename,'Tcline','units','meter')
% %--------------------------------------------------------------------------
nccreate(inifilename,'hc','Dimensions',{'one',dim_one},'FillValue','disable');
ncwriteatt(inifilename,'hc','long_name','S-coordinate parameter, critical depth')
ncwriteatt(inifilename,'hc','units','meter')
% %--------------------------------------------------------------------------
nccreate(inifilename,'s_rho','Dimensions',{'s_rho',dim_s_rho},'FillValue','disable');
ncwriteatt(inifilename,'s_rho','long_name','S-coordinate at RHO-points')
ncwriteatt(inifilename,'s_rho','valid_min',-1)
ncwriteatt(inifilename,'s_rho','valid_max',0)
ncwriteatt(inifilename,'s_rho','positive','up')
ncwriteatt(inifilename,'s_rho','standard_name','ocean_s_coordinate_g1')
ncwriteatt(inifilename,'s_rho','formula_terms','s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc')
% %--------------------------------------------------------------------------
nccreate(inifilename,'s_w','Dimensions',{'s_w',dim_s_w},'FillValue','disable');
ncwriteatt(inifilename,'s_w','long_name','S-coordinate at W-points')
ncwriteatt(inifilename,'s_w','valid_min',-1)
ncwriteatt(inifilename,'s_w','valid_max',0)
ncwriteatt(inifilename,'s_w','positive','up')
ncwriteatt(inifilename,'s_w','standard_name','ocean_s_coordinate_g1')
ncwriteatt(inifilename,'s_w','formula_terms','s: s_w C: Cs_w eta: zeta depth: h depth_c: hc')
% %--------------------------------------------------------------------------
nccreate(inifilename,'Cs_r','Dimensions',{'s_rho',dim_s_rho},'FillValue','disable');
ncwriteatt(inifilename,'Cs_r','long_name','S-coordinate stretching function at RHO-points')
ncwriteatt(inifilename,'Cs_r','valid_min',-1)
ncwriteatt(inifilename,'Cs_r','valid_max',0)
% %--------------------------------------------------------------------------
nccreate(inifilename,'Cs_w','Dimensions',{'s_w',dim_s_w},'FillValue','disable');
ncwriteatt(inifilename,'Cs_w','long_name','S-coordinate stretching function at W-points')
ncwriteatt(inifilename,'Cs_w','valid_min',-1)
ncwriteatt(inifilename,'Cs_w','valid_max',0)
% %--------------------------------------------------------------------------
nccreate(inifilename,'h','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho},'FillValue','disable');
ncwriteatt(inifilename,'h','long_name','bathymetry at RHO-points')
ncwriteatt(inifilename,'h','units','meter')
% %--------------------------------------------------------------------------
nccreate(inifilename,'lon_rho','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho},'FillValue','disable');
ncwriteatt(inifilename,'lon_rho','long_name','longitude of RHO-points')
ncwriteatt(inifilename,'lon_rho','units','degree_east')
% %--------------------------------------------------------------------------
nccreate(inifilename,'lat_rho','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho},'FillValue','disable');
ncwriteatt(inifilename,'lat_rho','long_name','latitute of RHO-points')
ncwriteatt(inifilename,'lat_rho','units','degree_north')
% %--------------------------------------------------------------------------
nccreate(inifilename,'lon_u','Dimensions',{'xi_u',dim_xi_u,'eta_u',dim_eta_u},'FillValue','disable');
ncwriteatt(inifilename,'lon_u','long_name','longitude of U-points')
ncwriteatt(inifilename,'lon_u','units','degree_east')
% %--------------------------------------------------------------------------
nccreate(inifilename,'lat_u','Dimensions',{'xi_u',dim_xi_u,'eta_u',dim_eta_u},'FillValue','disable');
ncwriteatt(inifilename,'lat_u','long_name','latitute of U-points')
ncwriteatt(inifilename,'lat_u','units','degree_north')
% %--------------------------------------------------------------------------
nccreate(inifilename,'lon_v','Dimensions',{'xi_v',dim_xi_v,'eta_v',dim_eta_v},'FillValue','disable');
ncwriteatt(inifilename,'lon_v','long_name','longitude of V-points')
ncwriteatt(inifilename,'lon_v','units','degree_east')
% %--------------------------------------------------------------------------
nccreate(inifilename,'lat_v','Dimensions',{'xi_v',dim_xi_v,'eta_v',dim_eta_v},'FillValue','disable');
ncwriteatt(inifilename,'lat_v','long_name','latitute of V-points')
ncwriteatt(inifilename,'lat_v','units','degree_north')
%--------------------------------------------------------------------------
nccreate(inifilename,'ocean_time','Dimensions',{'ocean_time',dim_ocean_time},'FillValue','disable');
ncwriteatt(inifilename,'ocean_time','long_name','time since initialization')
ncwriteatt(inifilename,'ocean_time','units','second since 2003-01-01 00:00:00')
ncwriteatt(inifilename,'ocean_time','calendar','365.0 days in every year')
% %--------------------------------------------------------------------------
nccreate(inifilename,'zeta','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
ncwriteatt(inifilename,'zeta','long_name','free-surface')
ncwriteatt(inifilename,'zeta','units','meter')
ncwriteatt(inifilename,'zeta','time','ocean_time')
ncwriteatt(inifilename,'zeta','coordinates','lon_rho lat_rho ocean_time')
% %--------------------------------------------------------------------------
nccreate(inifilename,'ubar','Dimensions',{'xi_u',dim_xi_u,'eta_u',dim_eta_u,'ocean_time',dim_ocean_time},'FillValue','disable');
ncwriteatt(inifilename,'ubar','long_name','vertically integrated u-momentum component')
ncwriteatt(inifilename,'ubar','units','meter second-1')
ncwriteatt(inifilename,'ubar','time','ocean_time')
ncwriteatt(inifilename,'ubar','coordinates','lon_u lat_u ocean_time')
% %--------------------------------------------------------------------------
nccreate(inifilename,'vbar','Dimensions',{'xi_v',dim_xi_v,'eta_v',dim_eta_v,'ocean_time',dim_ocean_time},'FillValue','disable');
ncwriteatt(inifilename,'vbar','long_name','vertically integrated v-momentum component')
ncwriteatt(inifilename,'vbar','units','meter second-1')
ncwriteatt(inifilename,'vbar','time','ocean_time')
ncwriteatt(inifilename,'vbar','coordinates','lon_v lat_v ocean_time')
% %--------------------------------------------------------------------------
nccreate(inifilename,'u','Dimensions',{'xi_u',dim_xi_u,'eta_u',dim_eta_u,'s_rho',dim_s_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
ncwriteatt(inifilename,'u','long_name','u-momentum component')
ncwriteatt(inifilename,'u','units','meter second-1')
ncwriteatt(inifilename,'u','time','ocean_time')
ncwriteatt(inifilename,'u','coordinates','lon_u lat_u s_rho ocean_time')
% %--------------------------------------------------------------------------
nccreate(inifilename,'v','Dimensions',{'xi_v',dim_xi_v,'eta_v',dim_eta_v,'s_rho',dim_s_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
ncwriteatt(inifilename,'v','long_name','v-momentum component')
ncwriteatt(inifilename,'v','units','meter second-1')
ncwriteatt(inifilename,'v','time','ocean_time')
ncwriteatt(inifilename,'v','coordinates','lon_v lat_v s_rho ocean_time')
% %--------------------------------------------------------------------------
nccreate(inifilename,'temp','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'s_rho',dim_s_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
ncwriteatt(inifilename,'temp','long_name','potential temperature')
ncwriteatt(inifilename,'temp','units','Celsius')
ncwriteatt(inifilename,'temp','time','ocean_time')
ncwriteatt(inifilename,'temp','coordinates','lon_rho lat_rho s_rho ocean_time')
% %--------------------------------------------------------------------------
nccreate(inifilename,'salt','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'s_rho',dim_s_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
ncwriteatt(inifilename,'salt','long_name','salinity')
ncwriteatt(inifilename,'salt','units','PSU')
ncwriteatt(inifilename,'salt','time','ocean_time')
ncwriteatt(inifilename,'salt','coordinates','lon_rho lat_rho s_rho ocean_time')
% %--------------------------------------------------------------------------
nccreate(inifilename,'ageice','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'uice','Dimensions',{'xi_u',dim_xi_u,'eta_u',dim_eta_u,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'vice','Dimensions',{'xi_v',dim_xi_v,'eta_v',dim_eta_v,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'aice','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'hice','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'tisrf','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'snow_thick','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'sfwat','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'ti','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'sig11','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'sig22','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'sig12','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'tau_iw','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'chu_iw','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'s0mk','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'t0mk','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'dye_01','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'s_rho',dim_s_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
% %--------------------------------------------------------------------------
nccreate(inifilename,'dye_02','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'s_rho',dim_s_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
%--------------------------------------------------------------------------
nccreate(inifilename,'dye_03','Dimensions',{'xi_rho',dim_xi_rho,'eta_rho',dim_eta_rho,'s_rho',dim_s_rho,'ocean_time',dim_ocean_time},'FillValue','disable');
%--------------------------------------------------------------------------
% --- fill Variables ---
disp('    => Filling Variables...')
% 
%--------------------------------------------------------------------------
inivarfilename='F:\Ross_amundsen_roms_model\initialization_file\data\RAS_ini.mat';
load(inivarfilename);
zeta_ini = zeros(Lp,Mp);
ubar_ini = zeros(L,Mp);
vbar_ini = zeros(Lp,M);
u_ini = zeros(L,Mp,N);
v_ini = zeros(Lp,M,N);
% temp_ini = zeros(Lp,Mp,N);
% salt_ini = zeros(Lp,Mp,N);
ageice_ini = zeros(Lp,Mp);
uice_ini = zeros(L,Mp);
vice_ini = zeros(Lp,M);
aice_ini = zeros(Lp,Mp);
hice_ini = zeros(Lp,Mp);
tisrf_ini = zeros(Lp,Mp);
snow_thick_ini = zeros(Lp,Mp);
sfwat_ini = zeros(Lp,Mp);
ti_ini = zeros(Lp,Mp);
sig11_ini = zeros(Lp,Mp);
sig22_ini = zeros(Lp,Mp);
sig12_ini = zeros(Lp,Mp);
tau_iw_ini = zeros(Lp,Mp);
chu_iw_ini = zeros(Lp,Mp);
s0mk_ini = zeros(Lp,Mp);
t0mk_ini = zeros(Lp,Mp);
% dye_01_ini = zeros(Lp,Mp,N);
dye_02_ini = zeros(Lp,Mp,N);
dye_03_ini = zeros(Lp,Mp,N);
%--------------------------------------------------------------------------
ncwrite(inifilename,'spherical',spherical);
ncwrite(inifilename,'Vtransform',Vtransform);
ncwrite(inifilename,'Vstretching',Vstretching);
ncwrite(inifilename,'theta_s',theta_s);
ncwrite(inifilename,'theta_b',theta_b);
ncwrite(inifilename,'Tcline',Tcline);
ncwrite(inifilename,'hc',hc);
ncwrite(inifilename,'s_rho',s_rho);
ncwrite(inifilename,'s_w',s_w);
ncwrite(inifilename,'Cs_r',Cs_r);
ncwrite(inifilename,'Cs_w',Cs_w);
ncwrite(inifilename,'h',h);
ncwrite(inifilename,'lon_rho',lon_rho);
ncwrite(inifilename,'lat_rho',lat_rho);
ncwrite(inifilename,'lon_u',lon_u);
ncwrite(inifilename,'lat_u',lat_u);
ncwrite(inifilename,'lon_v',lon_v);
ncwrite(inifilename,'lat_v',lat_v);
ncwrite(inifilename,'ocean_time',ocean_time_ini);
ncwrite(inifilename,'zeta',zeta_ini);
ncwrite(inifilename,'ubar',ubar_ini);
ncwrite(inifilename,'vbar',vbar_ini);
ncwrite(inifilename,'u',u_ini);
ncwrite(inifilename,'v',v_ini);
ncwrite(inifilename,'temp',temp_ini);
ncwrite(inifilename,'salt',salt_ini);
ncwrite(inifilename,'ageice',ageice_ini);
ncwrite(inifilename,'uice',uice_ini);
ncwrite(inifilename,'vice',vice_ini);
ncwrite(inifilename,'aice',aice_ini);
ncwrite(inifilename,'hice',hice_ini);
ncwrite(inifilename,'tisrf',tisrf_ini);
ncwrite(inifilename,'snow_thick',snow_thick_ini);
ncwrite(inifilename,'sfwat',sfwat_ini);
ncwrite(inifilename,'ti',ti_ini);
ncwrite(inifilename,'sig11',sig11_ini);
ncwrite(inifilename,'sig22',sig22_ini);
ncwrite(inifilename,'sig12',sig12_ini);
ncwrite(inifilename,'tau_iw',tau_iw_ini);
ncwrite(inifilename,'chu_iw',chu_iw_ini);
ncwrite(inifilename,'s0mk',s0mk_ini);
ncwrite(inifilename,'t0mk',t0mk_ini);
ncwrite(inifilename,'dye_01',dye_01_ini);
ncwrite(inifilename,'dye_02',dye_02_ini);
ncwrite(inifilename,'dye_03',dye_03_ini);
% ------------ write Global Attributes ----------
NLM_LBC(1,:)='EDGE   WEST   SOUTH  EAST   NORTH   ';
NLM_LBC(2,:)='zeta:  Cha    Clo    Cha    Cha     ';
NLM_LBC(3,:)='ubar:  Fla    Clo    Fla    Fla     ';
NLM_LBC(4,:)='vbar:  Fla    Clo    Fla    Fla     ';
NLM_LBC(5,:)='u:     Rad    Clo    Rad    Rad     ';
NLM_LBC(6,:)='v:     Rad    Clo    Rad    Rad     ';
NLM_LBC(7,:)='temp:  RadNud Clo    RadNud RadNud  ';
NLM_LBC(8,:)='salt:  RadNud Clo    RadNud RadNud  ';
ncwriteatt(inifilename,'/','type','INITIALIZATION file for RAS model');
ncwriteatt(inifilename,'/','creation_date',datestr(now));
ncwriteatt(inifilename,'/','history','..\RAS_model\RAS_ini.nc');
ncwriteatt(inifilename,'/','NLM_LBC',NLM_LBC');
disp('  Done !');



