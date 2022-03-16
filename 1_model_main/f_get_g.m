function gg = f_get_g(input_r1,input_r2)

rr = [input_r1,input_r2];
rr = 1.3*(rr+1.4e-10);
r1 = rr(1);
r2 = rr(2);

epsi = 0.51e-9;    
Rh = 42e-9;
aa = 0.53;
if (r1>r2)
  rp = r2;
else
  rp = r1;
end
gg = ((epsi/Rh)^2 + (epsi/rp)^2)^(-0.5*aa) / 0.22;
