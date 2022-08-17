%Process NOAANSIDC CDR of Passive Microwave sea ice concentration data
% RAS模型边界场文件中海冰密集度插值与处理
%Created by xiechuan  2021/10/20
clc;
clear all;
close all;
%
year_begin=1998;
year_end=1998;
%
ROMS_grdfiles_dir=['G:\Ross_amundsen_roms_model\grid_file\'];
grdfilename=[ROMS_grdfiles_dir,'RAS_grd_32layer_new.nc'];     
mask_rho=ncread(grdfilename,'mask_rho');
lon_rho=ncread(grdfilename,'lon_rho');
lat_rho=ncread(grdfilename,'lat_rho');
for i=1:1000
    for j=1:1000
        if lon_rho(i,j)<0
            lon_rho(i,j)=lon_rho(i,j)+360;
        end
    end
end
lon_west=lon_rho(1,:);
lat_west=lat_rho(1,:);
mask_west=mask_rho(1,:);
lon_east=lon_rho(end,:);
lat_east=lat_rho(end,:);
mask_east=mask_rho(end,:);
lon_north=lon_rho(:,end);
lat_north=lat_rho(:,end);
mask_north=mask_rho(:,end);

AMSR2_filedir='G:\Ross_amundsen_roms_model\boundary_file\A2_AMSR_sic_rawdata\AMSR2\';
AMSRE_filedir='G:\Ross_amundsen_roms_model\boundary_file\A2_AMSR_sic_rawdata\AMSRE\';
SSMIS_filedir='G:\Ross_amundsen_roms_model\boundary_file\A2_AMSR_sic_rawdata\SSMIS\';
SMMR_filedir='G:\Ross_amundsen_roms_model\boundary_file\A3_SMMR_sic_rawdata\';
output_datafir='G:\Ross_amundsen_roms_model\boundary_file\B2_Glosea5_after_interp\';

for iyear=year_begin:year_end
    %------------1998 to 2002-----------
    if iyear>=1998&&iyear<=2002
       % NOAA/NSIDC CDR data
       % NSIDC-processed Bootstrap daily sea ice concentrations      
        filename=[SMMR_filedir,'year\seaice_conc_daily_sh_1998_v04r00.nc'];
        lat_smmr=ncread(filename,'latitude');
        lon_smmr=ncread(filename,'longitude');
        lon_smmr(find(lon_smmr<0))=lon_smmr(find(lon_smmr<0))+360;
        for jmonth=1:12
           total_day=1;
           for kday=1:eomday(iyear,jmonth)
               filename=[SMMR_filedir,num2str(iyear),'\seaice_conc_daily_sh_',num2str(iyear),num2str(jmonth,'%02d'),num2str(kday,'%02d'),'_f13_v04r00.nc'];
               sic_tmp=ncread(filename,'nsidc_bt_seaice_conc');  
               sic_tmp(find(sic_tmp>1))=NaN;
               aice_east_tmp=griddata(lat_smmr,lon_smmr,sic_tmp,lat_east,lon_east,'cubic')';
               aice_west_tmp=griddata(lat_smmr,lon_smmr,sic_tmp,lat_west,lon_west,'cubic')';
               aice_north_tmp=griddata(lat_smmr,lon_smmr,sic_tmp,lat_north,lon_north,'cubic');
               aice_east_tmp(isnan(aice_east_tmp))=0;
               aice_west_tmp(isnan(aice_west_tmp))=0;
               aice_north_tmp(isnan(aice_north_tmp))=0;
               aice_east_tmp(find(aice_east_tmp<0))=0;
               aice_west_tmp(find(aice_west_tmp<0))=0;
               aice_north_tmp(find(aice_north_tmp<0))=0;
               aice_east(:,total_day)=aice_east_tmp.*mask_east';
               aice_west(:,total_day)=aice_west_tmp.*mask_west';
               aice_north(:,total_day)=aice_north_tmp.*mask_north;
               total_day=total_day+1;
           end
           bryfilename=[output_datafir,'RAS_bry_data_',num2str(iyear),num2str(jmonth,'%02d'),'.mat'];
           save(bryfilename,'aice_east','aice_west','aice_north','-append');
           clear aice_east aice_west aice_north
           disp(['Processing : ',num2str(iyear),num2str(jmonth,'%02d'),' ... ']);
        end
    end
    %----------------2003 to 2010--------------
    if iyear>=2003&&data_year<=2010
        % 处理每年的海冰密集度数据AMSRE       
        filename_grid=[AMSRE_filedir,'LongitudeLatitudeGrid-s3125-Antarctic3125.hdf'];
        lon_amsr=hdfread(filename_grid,'Longitudes');
        lat_amsr=hdfread(filename_grid,'Latitudes');
        for jmonth=1:12
           total_day=1;
           for kday=1:eomday(iyear,jmonth)
               filename=[AMSRE_filedir,num2str(iyear),'\asi-AMSRE-s3125-',num2str(iyear),num2str(jmonth,'%02d'),num2str(kday,'%02d'),'-v5.4.hdf'];
               sic_tmp=hdfread(filename,'ASI Ice Concentration');
               aice_east(:,total_day)=griddata(double(lat_amsr),double(lon_amsr),double(sic_tmp),lat_east,lon_east,'cubic')';
               aice_west(:,total_day)=griddata(double(lat_amsr),double(lon_amsr),double(sic_tmp),lat_west,lon_west,'cubic')';
               aice_north(:,total_day)=griddata(double(lat_amsr),double(lon_amsr),double(sic_tmp),lat_north,lon_north,'cubic');
               aice_east(isnan(aice_east(:,total_day)),total_day)=0;
               aice_west(isnan(aice_west(:,total_day)),total_day)=0;
               aice_north(isnan(aice_north(:,total_day)),total_day)=0;
               aice_east(:,total_day)=aice_east(:,total_day).*mask_east';
               aice_west(:,total_day)=aice_west(:,total_day).*mask_west';
               aice_north(:,total_day)=aice_north(:,total_day).*mask_north;
               total_day=total_day+1;
           end
           bryfilename=[output_datafir,'RAS_bry_data_',num2str(iyear),num2str(jmonth,'%02d'),'.mat'];
           save(bryfilename,'aice_east','aice_west','aice_north','-append');
           clear aice_east aice_west aice_north
           disp(['Processing : ',num2str(iyear),num2str(jmonth,'%02d'),' ... ']);
        end
    end
    %----------------------2011 to 2012-------------------
    if iyear>=2011&&iyear<=2012
        % 处理每年的海冰密集度数据SSMIS
        filename_grid=[SSMIS_filedir,'LongitudeLatitudeGrid-s6250-Antarctic.hdf'];
        lon_amsr=hdfread(filename_grid,'Longitudes');
        lat_amsr=hdfread(filename_grid,'Latitudes');
        for jmonth=1:12
           total_day=1;
           for kday=1:eomday(iyear,jmonth)
               filename=[AMSR2_filedir,num2str(iyear),'\asi-SSMIS17-s6250-',num2str(iyear),num2str(jmonth,'%02d'),num2str(kday,'%02d'),'-v5.hdf'];
               sic_tmp=hdfread(filename,'ASI Ice Concentration');
               aice_east(:,total_day)=griddata(double(lat_amsr),double(lon_amsr),double(sic_tmp),lat_east,lon_east,'cubic')';
               aice_west(:,total_day)=griddata(double(lat_amsr),double(lon_amsr),double(sic_tmp),lat_west,lon_west,'cubic')';
               aice_north(:,total_day)=griddata(double(lat_amsr),double(lon_amsr),double(sic_tmp),lat_north,lon_north,'cubic');
               aice_east(isnan(aice_east(:,total_day)),total_day)=0;
               aice_west(isnan(aice_west(:,total_day)),total_day)=0;
               aice_north(isnan(aice_north(:,total_day)),total_day)=0;
               aice_east(:,total_day)=aice_east(:,total_day).*mask_east';
               aice_west(:,total_day)=aice_west(:,total_day).*mask_west';
               aice_north(:,total_day)=aice_north(:,total_day).*mask_north;
               total_day=total_day+1;
           end
           bryfilename=[output_datafir,'RAS_bry_data_',num2str(iyear),num2str(jmonth,'%02d'),'.mat'];
           save(bryfilename,'aice_east','aice_west','aice_north','-append');
           clear aice_east aice_west aice_north
           disp(['Processing : ',num2str(iyear),num2str(jmonth,'%02d'),' ... ']);
        end
    end
    %--------------------2013 to 2016-------------------
    if iyear>=2013
        % 处理每年的海冰密集度数据AMSR2        
        filename_grid=[AMSR2_filedir,'LongitudeLatitudeGrid-s3125-Antarctic3125.hdf'];
        lon_amsr=hdfread(filename_grid,'Longitudes');
        lat_amsr=hdfread(filename_grid,'Latitudes');
        for jmonth=1:12
           total_day=1; 
           for kday=1:eomday(iyear,jmonth)
               filename=[AMSR2_filedir,num2str(iyear),'\asi-AMSR2-s3125-',num2str(iyear),num2str(jmonth,'%02d'),num2str(kday,'%02d'),'-v5.4.hdf'];
               sic_tmp=hdfread(filename,'ASI Ice Concentration');
               aice_east(:,total_day)=griddata(double(lat_amsr),double(lon_amsr),double(sic_tmp),lat_east,lon_east,'cubic')';
               aice_west(:,total_day)=griddata(double(lat_amsr),double(lon_amsr),double(sic_tmp),lat_west,lon_west,'cubic')';
               aice_north(:,total_day)=griddata(double(lat_amsr),double(lon_amsr),double(sic_tmp),lat_north,lon_north,'cubic');
               aice_east(isnan(aice_east(:,total_day)),total_day)=0;
               aice_west(isnan(aice_west(:,total_day)),total_day)=0;
               aice_north(isnan(aice_north(:,total_day)),total_day)=0;
               aice_east(:,total_day)=aice_east(:,total_day).*mask_east';
               aice_west(:,total_day)=aice_west(:,total_day).*mask_west';
               aice_north(:,total_day)=aice_north(:,total_day).*mask_north;
               total_day=total_day+1;
           end
           bryfilename=[output_datafir,'RAS_bry_data_',num2str(iyear),num2str(jmonth,'%02d'),'.mat'];
           save(bryfilename,'aice_east','aice_west','aice_north','-append');
           clear aice_east aice_west aice_north
           disp(['Processing : ',num2str(iyear),num2str(jmonth,'%02d'),' ... ']);
        end
    end
end




