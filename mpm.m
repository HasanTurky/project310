function cp = mpm(length, tx_h,rx_h, sf_h, order_n, r)
L = length;  %length in m
a = tx_h ; %transmitter height
b = rx_h ; %receiver height
h = sf_h; %surface height
o = order_n;
 %r surface reflection coefficient

D = L + 1/(2*L)*(b-a).^2;
SSn = L + 1/(2*L)*(2*o*h-b-a).^2;
SBn = L + 1/(2*L)*(2*o*h+b-a).^2;
BSn = L + 1/(2*L)*(2*o*h-b+a).^2;
BBn = L + 1/(2*L)*(2*(o-1)*h+a+b).^2;

%% combined attenuation
Rssn = -abs(r).^o;
Rsbn = abs(r).^o;
Rbsn = abs(r).^o;
Rbbn = -abs(r).^(o-1);

%% attenuation coefficient

a_ssn = D/SSn*Rssn;
a_sbn = D/SBn*Rsbn;
a_bsn = D/BSn*Rbsn;
a_bbn = D/BBn*Rbbn;

cp = [a_ssn a_sbn a_bsn a_bbn];
end
