function [cur_1overKM,complex_D,log_Gamma] = f_calc_d_without_ts(rho_vol,ttheta,arr_HSV_SED,vol_occ,ratio_s,N_protein,KM0)

[arr_H,arr_S,arr_V] = arr_HSV_SED{1:3};

xs0 = vol_occ*ratio_s/arr_V(1) / N_protein;
xE0 = vol_occ*(1-ratio_s)/arr_V(2) / N_protein;

%disp([xs0,xE0]);


%ffactor_radius = 1/arr_H(1)+1/arr_H(2);


xH = [arr_H;arr_H];
xS = [arr_S;arr_S];
xV = [arr_V;arr_V];



cur_1overKM = 1/KM0;
while 1
  cur_D = xE0/(1+1/cur_1overKM/xs0);
%disp([xE0,xs0,cur_1overKM,cur_D])

  
  x_tmp = [xs0-cur_D; xE0-cur_D; cur_D];
  x_count = [x_tmp; (N_protein-1)*x_tmp];
%disp(x_count');
  
  avg_V = sum(x_count.*xV)/rho_vol;
  avg_S = sum(x_count.*xS)/rho_vol;
  avg_H = sum(x_count.*xH)/rho_vol;
  avg_H2 = sum(x_count.*xH.*xH)/rho_vol;
  avg_1 = sum(x_count)/rho_vol;  
%disp([avg_V,avg_S,avg_H,avg_H2,avg_1]);
%disp([avg_V]);
  log_gamma = -1*log(1-avg_V) + (xH*avg_S + xS*avg_H + xV*avg_1) / (1-avg_V) + ((xH.^2)*avg_S^2 + 2*xV*avg_H*avg_S) / 2 / (1-avg_V)^2 + xV*avg_H2*avg_S^2 / 3 / (1-avg_V)^3;
%disp(log_gamma);
  
  %log_Gamma = log_gamma(1) + log_gamma(2) - log_gamma(3);
  log_Gamma = 1;
  
  pre_1overKM = cur_1overKM;
  
  if log_Gamma>500
%disp([avg_V log_Gamma]);
    %termm_denominator = (1+ttheta*ffactor_radius/const_ratio_radius)*exp(-1*g_diffusion*avg_V);
    termm_denominator = (1+ttheta)*exp(-1*f_get_g(arr_H(1),arr_H(2))*avg_V);
    termm_numerator = 1;
  else
    %termm_denominator = (1+ttheta*ffactor_radius/const_ratio_radius)*exp(log_Gamma-g_diffusion*avg_V);
    %termm_numerator = exp(log_Gamma) + ttheta*ffactor_radius/const_ratio_radius*exp(-1*g_diffusion*avg_V);
    termm_denominator = (1+ttheta)*exp(log_Gamma-f_get_g(arr_H(1),arr_H(2))*avg_V);
    termm_numerator = exp(log_Gamma) + ttheta*exp(-1*f_get_g(arr_H(1),arr_H(2))*avg_V);
  end
  cur_1overKM = 1/KM0 * termm_denominator / termm_numerator;
  
  complex_D = cur_D;
  if and(pre_1overKM==0,cur_1overKM==0)
    break
  elseif cur_1overKM==0
    ;
  else
    ttmp = abs(cur_1overKM-pre_1overKM)/pre_1overKM;
    if ttmp<1e-3
      break
    end
  end
end

