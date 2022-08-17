% 
clc;
clear all;

Year_begin=1998;
Year_end=1998;
bry_filename_path='G:\Ross_amundsen_roms_model\boundary_file\B2_Glosea5_after_interp\';
outputdir='G:\Ross_amundsen_roms_model\boundary_file\B3_RAS_boundary_varibles_yearly\';
for iyear=Year_begin:Year_end   
    bry_filename=[bry_filename_path,'RAS_bry_data_',num2str(iyear),'01.mat'];
    load(bry_filename);
    zeta_west_bry=zeta_west;
    zeta_east_bry=zeta_east;
    zeta_north_bry=zeta_north;
    salt_west_bry=salt_west;
    salt_east_bry=salt_east;
    salt_north_bry=salt_north;
    temp_west_bry=temp_west;
    temp_east_bry=temp_east;
    temp_north_bry=temp_north;
    ua_west_bry=ua_west;
    ua_east_bry=ua_east;
    ua_north_bry=ua_north;
    va_west_bry=va_west;
    va_east_bry=va_east;
    va_north_bry=va_north;
    aice_west_bry=aice_west;
    aice_east_bry=aice_east;
    aice_north_bry=aice_north;
    dye_west_01_bry=dye_west_01;
    dye_east_01_bry=dye_east_01;
    dye_north_01_bry=dye_north_01;   
    t_salt_bry=t_salt';
    t_zeta_bry=t_zeta';
    t_v2d_bry=t_v2d';
    t_temp_bry=t_temp';
    ts_bry=ts;
    for jmonth=2:12
        bry_filename=[bry_filename_path,'RAS_bry_data_',num2str(iyear),num2str(jmonth,'%02d'),'.mat'];
        load(bry_filename);
        zeta_west_bry=cat(2,zeta_west_bry,zeta_west);
        zeta_east_bry=cat(2,zeta_east_bry,zeta_east);
        zeta_north_bry=cat(2,zeta_north_bry,zeta_north);
        salt_west_bry=cat(3,salt_west_bry,salt_west);
        salt_east_bry=cat(3,salt_east_bry,salt_east);
        salt_north_bry=cat(3,salt_north_bry,salt_north);
        temp_west_bry=cat(3,temp_west_bry,temp_west);
        temp_east_bry=cat(3,temp_east_bry,temp_east);
        temp_north_bry=cat(3,temp_north_bry,temp_north);
%         u_west=cat(3,u_west,u_west);
%         u_east=cat(3,u_east,u_east);
%         u_north=cat(3,u_north,u_north);
%         v_west=cat(3,v_west,v_west);
%         v_east=cat(3,v_east,v_east);
%         v_north=cat(3,v_north,v_north);
        ua_west_bry=cat(2,ua_west_bry,ua_west);
        ua_east_bry=cat(2,ua_east_bry,ua_east);
        ua_north_bry=cat(2,ua_north_bry,ua_north);
        va_west_bry=cat(2,va_west_bry,va_west);
        va_east_bry=cat(2,va_east_bry,va_east);
        va_north_bry=cat(2,va_north_bry,va_north);
        aice_west_bry=cat(2,aice_west_bry,aice_west);
        aice_east_bry=cat(2,aice_east_bry,aice_east);
        aice_north_bry=cat(2,aice_north_bry,aice_north);
        dye_west_01_bry=cat(3,dye_west_01_bry,dye_west_01);
        dye_east_01_bry=cat(3,dye_east_01_bry,dye_east_01);
        dye_north_01_bry=cat(3,dye_north_01_bry,dye_north_01);
        t_salt_bry=cat(2,t_salt_bry,t_salt');
        t_zeta_bry=cat(2,t_zeta_bry,t_zeta');
        t_v2d_bry=cat(2,t_v2d_bry,t_v2d');
        t_temp_bry=cat(2,t_temp_bry,t_temp');
        ts_bry=cat(1,ts_bry,ts);
    end
    %----------add south boundary (zeros matrix)---------
    dye_south_01_bry=zeros(size(dye_north_01_bry));
    aice_south_bry=zeros(size(aice_north_bry));
    temp_south_bry=zeros(size(temp_north_bry));
    salt_south_bry=zeros(size(salt_north_bry));
    zeta_south_bry=zeros(size(zeta_north_bry));
    ua_south_bry=zeros(size(ua_north_bry));
    va_south_bry=zeros(size(va_north_bry));
    %---------add dye02&dye03 (zeros matrix)--------------
    dye_west_02_bry=zeros(size(dye_west_01_bry));
    dye_east_02_bry=zeros(size(dye_east_01_bry));
    dye_north_02_bry=zeros(size(dye_north_01_bry));
    dye_south_02_bry=zeros(size(dye_north_01_bry));
    dye_west_03_bry=zeros(size(dye_west_01_bry));
    dye_east_03_bry=zeros(size(dye_east_01_bry));
    dye_north_03_bry=zeros(size(dye_north_01_bry));
    dye_south_03_bry=zeros(size(dye_north_01_bry));
    %----------------save data-----------------
    filename=[outputdir,'RAS_bry_',num2str(iyear),'.mat'];
    save(filename,'zeta_west_bry','zeta_east_bry','zeta_north_bry','zeta_south_bry');
    save(filename,'salt_west_bry','salt_east_bry','salt_north_bry','salt_south_bry','-append');
    save(filename,'temp_west_bry','temp_east_bry','temp_north_bry','temp_south_bry','-append');
    save(filename,'ua_west_bry','ua_east_bry','ua_north_bry','ua_south_bry','-append');
    save(filename,'va_west_bry','va_east_bry','va_north_bry','va_south_bry','-append');
    save(filename,'dye_west_01_bry','dye_east_01_bry','dye_north_01_bry','dye_south_01_bry','-append');
    save(filename,'dye_west_02_bry','dye_east_02_bry','dye_north_02_bry','dye_south_02_bry','-append');
    save(filename,'dye_west_03_bry','dye_east_03_bry','dye_north_03_bry','dye_south_03_bry','-append');
    save(filename,'aice_west_bry','aice_east_bry','aice_north_bry','aice_south_bry','-append');
    save(filename,'t_salt_bry','t_zeta_bry','t_v2d_bry','t_temp_bry','ts_bry','-append');
    
end

