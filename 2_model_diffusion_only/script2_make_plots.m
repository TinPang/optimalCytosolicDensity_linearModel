arr_occupancy = 0.05:0.005:0.8;
arr_ratioS = 10.^(-3:0.001:-0.01);
mat_occ = arr_occupancy' * ones(1,length(arr_ratioS));
mat_ratioS = ones(length(arr_occupancy),1) * arr_ratioS;
mat_volS = mat_occ.*mat_ratioS;
mat_volR = mat_occ.*(1-mat_ratioS);
%contour(mat_volS,mat_volR,aa)



for Nprotein = 20:80:100

aa = load(['intermediateData_NComplex_smallEnzyme_KM078260_',int2str(Nprotein),'Proteins']);
bb = load(['intermediateData_NComplex_largeEnzyme_KM078260_',int2str(Nprotein),'Proteins']);
aa_K = load(['intermediateData_K_smallEnzyme_KM078260_',int2str(Nprotein),'Proteins']);
bb_K = load(['intermediateData_K_largeEnzyme_KM078260_',int2str(Nprotein),'Proteins']);


mat_volocc = mat_occ;
cca = aa./mat_volocc;
ccb = bb./mat_volocc;
%volocc = mat_volS+mat_volR;
volocc = mat_occ;
aaa = unique(sort(reshape(volocc,[],1)));
bbb = [];
bbb_jj = [];
for qq = 1:length(aaa)
  [ii,jj] = find(volocc==aaa(qq));
  tmpcc = zeros(length(ii),4);
  for qqq = 1:length(ii)
    tmpcc(qqq,1) = cca(ii(qqq),jj(qqq));
    tmpcc(qqq,2) = ccb(ii(qqq),jj(qqq));
    tmpcc(qqq,3) = aa_K(ii(qqq),jj(qqq));
    tmpcc(qqq,4) = bb_K(ii(qqq),jj(qqq));
  end
  tmp_bbb = zeros(1,4);
  tmp_jjj = zeros(1,2);  
  for qqq = 1:2
    [vv,kk] = max(tmpcc(:,qqq));
    tmp_bbb(qqq) = vv;
    tmp_bbb(qqq+2) = tmpcc(kk,qqq+2);
    tmp_jjj(qqq) = jj(kk);
  end
  bbb(qq,1:4) = tmp_bbb;
  bbb_jj(qq,1:2) = tmp_jjj;
end



ccc = bbb;
%figure; plot(aaa(1:end),ccc(:,1)/max(ccc(:,1)),'b-');
%hold on; plot(aaa(1:end),ccc(:,2)/max(ccc(:,2)),'r-');
%xlabel('occupancy \rho');
%ylabel('growth rate \mu (arbitrary scale)');
KM0 = 130*602;
% plots the K_M^*/K_M^0 vs occupancy curve, i.e. Fig 7 in the manuscript
figure; semilogy(aaa(1:end),ccc(:,3)/KM0,'b-');
hold on; semilogy(aaa(1:end),ccc(:,4)/KM0,'r-');
xlabel('occupancy \rho')
ylabel('K_M^* / K_M^0')

end



