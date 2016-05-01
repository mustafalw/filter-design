syms s z;
w = linspace(0,2,2001);
jw = linspace(0,2*1i,2001);
num = 0.2017;
B = 0.5371;
Wo = 0.808;
p = zeros([1 4]);
p(1) = -0.1222 + 0.9698i;
p(2) = -0.2949 + 0.4017i;
p(3) = conj(p(2));
p(4) = conj(p(1));
H_jw = zeros([1 2001]);
H_jw = num + H_jw;
D_s = 1;
for i = 1:4
    H_jw = H_jw./(jw - p(i));
    D_s = D_s*(s - p(i));
end
D_s = sym2poly(D_s);
figure(1);
plot(w,abs(H_jw));
N_s = num*((B*s)^4);
D_s = 1;
w = linspace(0,5,5001);
jw = linspace(0,5*1i,5001);
H_jw = num*((B.*jw).^4);
for i = 1:4
    H_jw = H_jw./(-w.^2 - p(i)*B.*jw + Wo^2);
    D_s = D_s*(s^2 -p(i)*B*s + Wo^2);
end
D_s = sym2poly(D_s);
N_s = sym2poly(N_s);
figure(2);
plot(w, abs(H_jw));
w = linspace(0,pi,1001);
jw = linspace(0,pi*1i,1001);
e_jw = exp(jw);
H_e_jw = num*((B.*(e_jw.^2 - 1)).^4);
D_z = 1;
N_z = 0.2017*B^4*(z^2 - 1)^4;
for i = 1:4
    H_e_jw = H_e_jw./((e_jw.^2).*(1 - p(i)*B + Wo^2) + 2.*e_jw.*(Wo^2 - 1) + (1 + p(i)*B + Wo^2));
    D_z = D_z*((z^2)*(1 - p(i)*B + Wo^2) + 2*z*(Wo^2 - 1) + (1 + p(i)*B + Wo^2));
end
D_z = sym2poly(D_z);
N_z = sym2poly(N_z);
figure(3);
w = w./pi;
plot(w, abs(H_e_jw));