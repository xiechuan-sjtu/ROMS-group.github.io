%  First step 1: Read and organize the Glosea5 variables
%  Ross sea and Amundsen Sea Model (RAS)
%  Written by Xie Chuan
%  Create date: 2021-9-10
%
clc;
clear all; 
close all; 

%-------------------------------------------------------------------------------
%--------------------------- generate Glosea5 data -----------------------------
%-------------------------------------------------------------------------------
Year_begin=1999;
Year_end=2016;
% Common parameters
mybasedate = datenum(1998,1,1);

Monthabb=['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
for iyear=Year_begin:Year_end
   for jmonth=1:12
       filename_a=['G:\Ross_amundsen_roms_model\boundary_file\A1_Glosea5_rawdata\'...
           ,num2str(iyear),'\Glosea5_',num2str(iyear),num2str(jmonth,'%02d'),'_a.nc'];
       filename_b=['G:\Ross_amundsen_roms_model\boundary_file\A1_Glosea5_rawdata\'...
           ,num2str(iyear),'\Glosea5_',num2str(iyear),num2str(jmonth,'%02d'),'_b.nc'];
       salinity=ncread(filename_a,'so_foam');
       temp=ncread(filename_a,'thetao_foam');
       u=ncread(filename_b,'uo_foam');
       v=ncread(filename_b,'vo_foam');
       ssh=ncread(filename_b,'zos_foam');
       %Change the layer number to 1:72
       salinity=salinity(:,:,1:72,:);
       temp=temp(:,:,1:72,:);
       u=u(:,:,1:72,:);
       v=v(:,:,1:72,:);
       %date & day number
       num_days=eomday(iyear,jmonth); 
       for kdays=1:num_days
           ts(kdays,:)=[num2str(kdays,'%02d'),'-',Monthabb(jmonth,:),'-',num2str(iyear)];
           td(kdays,:)=datenum(ts(kdays,:))-mybasedate;
       end       
       filename_total=['G:\Ross_amundsen_roms_model\boundary_file\B1_Glosea5_before_interp\Glosea5_'...
           ,num2str(iyear),num2str(jmonth,'%02d'),'.mat'];
       save(filename_total,'salinity','temp','u','v','ssh','ts','td');
       clear ts td
       clear salinity temp u v ssh
       disp(['Processing : ',num2str(iyear),num2str(jmonth,'%02d'),' ... ']);
   end
end

%%
%更改网格
depth=depth(1:72,:);
save('Glosea5_grid_ross.mat','depth','Glosea5_mask','lat_Glosea5','lon_Glosea5');








