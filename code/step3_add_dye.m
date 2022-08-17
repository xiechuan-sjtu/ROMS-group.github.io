%% Add CDW dye (dye_01) in the RAS model bounary file
%  Written by Xie Chuan
%  Create date: 2021-9-10
%
%  Definition Rule （Dinniman and Klinck et al. 2012）
% (1) temperature greater than 0.0 dec below 200 m
% (2) placed off the continental shelf (no dye was initialized on the shelf (h < 1000m )) 
% (3) initial concentration of 1.

clc
clear all;
close all;

Year_begin = 1998;
Year_end = 1998;
grdfilename='G:\Ross_amundsen_roms_model\grid_file\RAS_grd_32layer_new.nc';
bry_filename_path='G:\Ross_amundsen_roms_model\boundary_file\B2_Glosea5_after_interp\';
load('G:\Ross_amundsen_roms_model\grid_file\mask_shelf.mat');
load('G:\Ross_amundsen_roms_model\grid_file\RAS_depth.mat');
mask_rho=ncread(grdfilename,'mask_rho');
% --- select open boundaries (by true/false) ---
obc_east=true;
obc_west=true;
obc_north=true;
obc_south=false;
%------------------------------------------------
for iyear=Year_begin:Year_end    
    for jmonth=1:12
        bry_filename=[bry_filename_path,'RAS_bry_data_',num2str(iyear),num2str(jmonth,'%02d'),'.mat'];
        load(bry_filename);
        size_east=size(temp_east);
        size_west=size(temp_west);
        size_north=size(temp_north);
        dye_east_01=zeros(size_east);
        dye_west_01=zeros(size_west);
        dye_north_01=zeros(size_north);
        zr_east=(-squeeze(zr(end,:,:)))';     %东边界
        zr_west=(-squeeze(zr(1,:,:)))';       %西边界
        zr_north=(-squeeze(zr(:,end,:)))';    %北边界
        %------------east-------------
        if obc_east
            for  i=1:size_east(3)
                for j=1:size_east(2)
                   for k=1:size_east(1)
                       if mask_rho(end,k)==1 && mask_shelf(end,k)==0
                           if temp_east(k,j,i)>0 && zr_east(j,k)>200
                              dye_east_01(k,j,i)=1;
                           end
                       end
                   end
                end
            end
            save(bry_filename,'dye_east_01','-append');
        end
        %--------------------west-----------------
        if obc_west
            for  i=1:size_west(3)
                for j=1:size_west(2)
                   for k=1:size_west(1)                  
                       if mask_rho(1,k)==1 && mask_shelf(1,k)==0
                           if temp_west(k,j,i)>0 && zr_west(j,k)>200
                               dye_west_01(k,j,i)=1;
                           end
                       end
                   end
                end
            end
            save(bry_filename,'dye_west_01','-append');
        end
        %------------------north----------------------
        if obc_north
            for  i=1:size_north(3)
                for j=1:size_north(2)
                       for k=1:size_north(1)
                           if mask_rho(k,end)==1 && mask_shelf(k,end)==0
                               if temp_north(k,j,i)>0 && zr_north(j,k)>200
                                   dye_north_01(k,j,i)=1;
                               end
                           end
                       end
                end
            end  
            save(bry_filename,'dye_north_01','-append');
        end
        disp(['Complete add CDW dye in ',num2str(iyear),num2str(jmonth,'%02d'),'...']);
    end
end
