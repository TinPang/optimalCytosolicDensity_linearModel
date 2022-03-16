%KM0 = 120400;
KM0 = 130*602;

E_r = 2.4e-9;
s_r = 0.34e-9;
%E_r = 10e-9;
%s_r = 2.4e-9;
rho_vol = 1e-18;

%g_diffusion = 6;
ttheta = 2.3;
%const_ratio_radius = 2/2e-9;



D_r = (E_r^3+s_r^3)^(1/3);
tmp_r = [s_r,E_r,D_r]';

arr_H = tmp_r;
arr_S = 4*pi*tmp_r.^2;
arr_V = 4/3*pi*tmp_r.^3;

%misc_r = [s_r,E_r,D_r]';
%misc_H = (misc_r*tan(pi/3) + misc_r/2*(pi*4/3)) / 4;
%misc_S = pi.*misc_r.^2 + pi.*misc_r.^2./cos(pi/3);
%misc_V = pi*misc_r.^3/3*tan(pi/3);
%misc_ratio = 1 ./ misc_V;
%misc_ratio = misc_ratio / sum(misc_ratio);

arr_HSV_SED = {arr_H,arr_S,arr_V};
%arr_HSV_misc = {misc_H,misc_S,misc_V,misc_ratio};
ratio_S2E = 1;

arr_occupancy = 0.05:0.005:0.8;
%arr_N_misc = rho_vol * arr_occupancy_misc / sum(misc_ratio.*misc_V);
arr_ratio_s = 10.^(-3:0.001:-0.01);



%%%%
%[cur_1overKM,complex_D,log_Gamma] = f_calc_c(rho_vol,g_diffusion,ttheta,const_ratio_radius,arr_HSV_SED,rho_vol*0.5,0.5,N_protein,KM0);


rec_occ_ratioS_maxD = [];
for N_protein = 20:80:100
arr_complex = zeros(length(arr_occupancy),length(arr_ratio_s));
arr_KM = zeros(length(arr_occupancy),length(arr_ratio_s));
arr_log_Gamma = zeros(length(arr_occupancy),length(arr_ratio_s));
for q = 1:length(arr_occupancy)
  %N_misc = arr_N_misc(q);
  for r = 1:length(arr_ratio_s)
      [cur_1overKM,complex_D,log_Gamma] = f_calc_d_without_G(rho_vol,ttheta,arr_HSV_SED,rho_vol*arr_occupancy(q),arr_ratio_s(r),N_protein,KM0);
    
      arr_complex(q,r) = complex_D;
      arr_KM(q,r) = 1/cur_1overKM;    
  end  
end
[mmax,imax_occ] = max(max(arr_complex,[],2));
[mmax,imax_ratioS] = max(max(arr_complex));
rec_occ_ratioS_maxD = [rec_occ_ratioS_maxD;[arr_occupancy(imax_occ),arr_ratio_s(imax_ratioS),arr_KM(imax_occ,imax_ratioS),arr_complex(imax_occ,imax_ratioS),mmax]];
dlmwrite(['tmp_ssave_NComplex_smallEnzyme_KM078260_',int2str(N_protein),'Proteins'],arr_complex);
dlmwrite(['tmp_ssave_K_smallEnzyme_KM078260_',int2str(N_protein),'Proteins'],arr_KM);
%dlmwrite(['tmp_ssave_NComplex_largeEnzyme_KM078260_',int2str(N_protein),'Proteins'],arr_complex);
%dlmwrite(['tmp_ssave_K_largeEnzyme_KM078260_',int2str(N_protein),'Proteins'],arr_KM);
end


%dlmwrite('tmp_ssave_NComplex_largeEnzyme_KM078260_20Proteins',arr_complex);
%dlmwrite('tmp_ssave_K_largeEnzyme_KM078260_20Proteins',arr_KM);

%dlmwrite('tmp_ssave_occ_ratioS_KM_maxD_5proteins_large',rec_occ_ratioS_maxD);



%20 proteins, (E_r,s_r) = (2e-9,0.43e-9), KM0 = 78260, optimal_occ = 55%




%arr_D = max(arr_complex,[],2);

%[cur_1overKM,complex_D,log_Gamma] = f_calc_c(rho_vol,g_diffusion,ttheta,const_ratio_radius,arr_HSV_SED,rho_vol*arr_occupancy(q),arr_ratio_s(r),N_protein,KM0);



% 20 proteins: 0.50
% 100 proteins: 0.40
