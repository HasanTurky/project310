function cp = mpm(length, tx_h,rx_h, sf_h, order_n, r)

D = L + 1/(2*L)*(b-a)^2;
SSn = L + 1/(2*L)*(2*n*h-b-a)^2;
SBn = L + 1/(2*L)*(2*n*h+b-a)^2;
BSn = L + 1/(2*L)*(2*n*h-b+a)^2;
BBn = L + 1/(2*L)*(2*(n-1)*h+a+b)^2;

%% combined attenuation
Rssn = -abs(r)^n;
Rsbn = abs(r)^n;
Rbsn = abs(r)^n;
Rbbn = -abs(r)^(n-1);

%% attenuation coefficient

a_ssn = D/SSn*Rssn;
a_sbn = D/RSn*Rsbn;
a_bsn = D/BSn*Rbsn;
a_bbn = D/BBn*Rbbn;

cp = [a_ssn a_sbn a_bsn a_bbn];
end
