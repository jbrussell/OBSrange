function [ data_all,misfit,lgd,symbol ] = agg_synthsurveys( filenames,Nfils )
%Aggregate synthetic survey results
misfit.agg = [];
misfit.agg.dtwt_all = [];
misfit.agg.dtwtcorr_all = [];
misfit.agg.twtcorr_all = [];
for ifil = 1:Nfils
    load(filenames{ifil});
    
    lat_drop = data(1).drop(1);
    lon_drop = data(1).drop(2);
    z_drop = data(1).drop(3)*-1000;
    lats_ship = data(1).survlats;
    lons_ship = data(1).survlons;
    olon = lon_drop;
    olat = lat_drop;
    [ x_ship, y_ship ] = lonlat2xy_nomap( olon, olat, lons_ship, lats_ship );
    data_all(ifil).x_ship = x_ship;
    data_all(ifil).y_ship = y_ship;
    data_all(ifil).v_surv_true = data(1).v_surv_true;
    data_all(ifil).vmag = sqrt(sum(data(1).v_surv_true.^2,2));
    data_all(ifil).survx = data(1).survx;
    data_all(ifil).survy = data(1).survy;
    data_all(ifil).survey = data(1).survey;
    data_all(ifil).radius = data(1).radius;
    
    
    misfit.xsta(ifil,1) = rms(data(1).misfit_xsta);
    misfit.xsta_std(ifil,1) = std(data(1).misfit_xsta);
    misfit.ysta(ifil,1) = rms(data(1).misfit_ysta);
    misfit.ysta_std(ifil,1) = std(data(1).misfit_ysta);
    misfit.zsta(ifil,1) = rms(data(1).misfit_zsta);
    misfit.zsta_std(ifil,1) = std(data(1).misfit_zsta);
    misfit.r_xy(ifil,1) = rms(data(1).misfit_r_xy);
    misfit.r_xy_std(ifil,1) = std(data(1).misfit_r_xy);
    misfit.r_xyz(ifil,1) = rms(data(1).misfit_r_xyz);
    misfit.r_xyz_std(ifil,1) = std(data(1).misfit_r_xyz);
    misfit.TAT(ifil,1) = rms(data(1).misfit_TAT);
    misfit.TAT_std(ifil,1) = std(data(1).misfit_TAT);
    misfit.Vw(ifil,1) = rms(data(1).misfit_Vw);
    misfit.Vw_std(ifil,1) = std(data(1).misfit_Vw);
    misfit.E_rms(ifil,1) = mean(data(1).E_rms);
    misfit.E_rms_std(ifil,1) = std(data(1).E_rms);
    misfit.v_ship_all(1:2,ifil) = mean(data(1).misfit_v_ship_all,2);
    misfit.v_ship_all_std(1:2,ifil) = std(data(1).misfit_v_ship_all,0,2);
    misfit.dtwtcorr_all(ifil,1) = rms(data(1).misfit_dtwtcorr_all);
    misfit.dtwtcorr_all_std(ifil,1) = std(data(1).misfit_dtwtcorr_all);
    misfit.dtwt_all(ifil,1) = rms(data(1).dtwt_all);
    misfit.dtwt_all_std(ifil,1) = std(data(1).dtwt_all);
    misfit.agg.dtwt_all = [misfit.agg.dtwt_all data(1).dtwt_all];
    misfit.agg.dtwtcorr_all = [misfit.agg.dtwtcorr_all data(1).misfit_dtwtcorr_all];
    misfit.agg.twtcorr_all = [misfit.agg.twtcorr_all vertcat(data.dtwtcorr)'];
    
    lgd{ifil} = [data_all(ifil).survey,' ',num2str(data_all(ifil).radius),' nm']; %files(ifil).name(1:end-4);
    if any(regexp(lgd{ifil},'PACMAN'))
        symbol{ifil} = 'ok';
    elseif any(regexp(lgd{ifil},'cross'))
        symbol{ifil} = 'sk';
    elseif any(regexp(lgd{ifil},'diamond'))
        symbol{ifil} = 'dk';
    elseif any(regexp(lgd{ifil},'line'))
        symbol{ifil} = 'pk';
    elseif any(regexp(lgd{ifil},'tri'))
        symbol{ifil} = '^k';
    end
    
end

end

