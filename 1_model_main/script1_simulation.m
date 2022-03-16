
% part 1: simulate small enzyme and small substrate (metabolic metabolic)

KM0 = 130*602;

E_r = 2.4e-9;  % enzyme radius; metric unit
s_r = 0.34e-9;  % substrate radius; metric unit
rho_vol = 1e-18;  % volume of a cell; metric unit

ttheta = 2.3;

D_r = (E_r^3+s_r^3)^(1/3);  % radius of enzyme-substrate-complex
tmp_r = [s_r,E_r,D_r]';
arr_H = tmp_r;  % stores r_i of Eq. S3
arr_S = 4*pi*tmp_r.^2;  % stores S_i of Eq. S3
arr_V = 4/3*pi*tmp_r.^3;  % stores V_i of Eq. S3

arr_HSV_SED = {arr_H,arr_S,arr_V};
ratio_S2E = 1;

arr_occupancy = 0.05:0.005:0.8;  % range of occupancy rho to consider
arr_ratio_s = 10.^(-3:0.001:-0.01);  % range of volume ratio (between enzyme and substrate) to consider

rec_occ_ratioS_maxD = [];
for N_protein = 20:80:100  % 'N_protein' corresponds to the number of protein 'N' in the linear model and parallel model; see Fig. 2
arr_complex = zeros(length(arr_occupancy),length(arr_ratio_s));
arr_KM = zeros(length(arr_occupancy),length(arr_ratio_s));
arr_log_Gamma = zeros(length(arr_occupancy),length(arr_ratio_s));
for q = 1:length(arr_occupancy)
  %N_misc = arr_N_misc(q);
  for r = 1:length(arr_ratio_s)
      [cur_1overKM,complex_D,log_Gamma] = f_calc_d_varyG(rho_vol,ttheta,arr_HSV_SED,rho_vol*arr_occupancy(q),arr_ratio_s(r),N_protein,KM0); % f_calc_d_varyG calculates the scaling factor g(r_s)    
      arr_complex(q,r) = complex_D;
      arr_KM(q,r) = 1/cur_1overKM;    
  end  
end
[mmax,imax_occ] = max(max(arr_complex,[],2));
[mmax,imax_ratioS] = max(max(arr_complex));
rec_occ_ratioS_maxD = [rec_occ_ratioS_maxD;[arr_occupancy(imax_occ),arr_ratio_s(imax_ratioS),arr_KM(imax_occ,imax_ratioS),arr_complex(imax_occ,imax_ratioS),mmax]];
dlmwrite(['intermediateData_NComplex_smallEnzyme_KM078260_',int2str(N_protein),'Proteins'],arr_complex);
dlmwrite(['intermediateData_K_smallEnzyme_KM078260_',int2str(N_protein),'Proteins'],arr_KM);
end




% part 2: simulate large ribosome and large tRNA (ribosomal system)

KM0 = 130*602;

E_r = 10e-9;  % enzyme radius; metric unit
s_r = 2.4e-9;  % substrate radius; metric unit
rho_vol = 1e-18;  % volume of a cell; metric unit

ttheta = 2.3;

D_r = (E_r^3+s_r^3)^(1/3);  % radius of enzyme-substrate-complex
tmp_r = [s_r,E_r,D_r]';
arr_H = tmp_r;  % stores r_i of Eq. S3
arr_S = 4*pi*tmp_r.^2;  % stores S_i of Eq. S3
arr_V = 4/3*pi*tmp_r.^3;  % stores V_i of Eq. S3

arr_HSV_SED = {arr_H,arr_S,arr_V};
ratio_S2E = 1;

arr_occupancy = 0.05:0.005:0.8;  % range of occupancy rho to consider
arr_ratio_s = 10.^(-3:0.001:-0.01);  % range of volume ratio (between enzyme and substrate) to consider

rec_occ_ratioS_maxD = [];
for N_protein = 20:80:100  % 'N_protein' corresponds to the number of protein 'N' in the linear model and parallel model; see Fig. 2
arr_complex = zeros(length(arr_occupancy),length(arr_ratio_s));
arr_KM = zeros(length(arr_occupancy),length(arr_ratio_s));
arr_log_Gamma = zeros(length(arr_occupancy),length(arr_ratio_s));
for q = 1:length(arr_occupancy)
  %N_misc = arr_N_misc(q);
  for r = 1:length(arr_ratio_s)
      [cur_1overKM,complex_D,log_Gamma] = f_calc_d_varyG(rho_vol,ttheta,arr_HSV_SED,rho_vol*arr_occupancy(q),arr_ratio_s(r),N_protein,KM0); % f_calc_d_varyG calculates the scaling factor g(r_s)    
      arr_complex(q,r) = complex_D;
      arr_KM(q,r) = 1/cur_1overKM;    
  end  
end
[mmax,imax_occ] = max(max(arr_complex,[],2));
[mmax,imax_ratioS] = max(max(arr_complex));
rec_occ_ratioS_maxD = [rec_occ_ratioS_maxD;[arr_occupancy(imax_occ),arr_ratio_s(imax_ratioS),arr_KM(imax_occ,imax_ratioS),arr_complex(imax_occ,imax_ratioS),mmax]];
dlmwrite(['intermediateData_NComplex_largeEnzyme_KM078260_',int2str(N_protein),'Proteins'],arr_complex);
dlmwrite(['intermediateData_K_largeEnzyme_KM078260_',int2str(N_protein),'Proteins'],arr_KM);
end
