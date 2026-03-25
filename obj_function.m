function y = obj_function(alpha,beta,th,phi,alpha_v,beta_v,th_v,phi_v,k1,k2,h,theta,m,n)

tao=0.5;
p1=1/(0.2/2.6*(2*pi)^(1/2))*exp(-(alpha-0.2)^2/(2*(0.2/2.6)^2));
p2=1/(0.4/2.6*(2*pi)^(1/2))*exp(-(beta-0.4)^2/(2*(0.4/2.6)^2));
p3=1/(0.25*(2*pi)^(1/2))*exp(-(th-1.5)^2/(2*(0.25)^2));
p4=1/(0.25*(2*pi)^(1/2))*exp(-(phi-1)^2/(2*(0.25)^2));

s1=0.001:0.05:20;
s=s1*1i;
num=(beta*s+alpha/th);
den=s.^2.*exp(phi*s)+(alpha/th+beta)*s+alpha/th;
H=num./den;
num=(beta_v*s+alpha_v/th_v);
den=s.^2.*exp(phi_v*s)+(alpha_v/th_v+beta_v)*s+alpha_v/th_v;
H_v=num./den;

D=tao*s+1;
T=1./(1+h*s);
G=D.*T.*(H_v).^m;
num=k2*s+k1+G.*s.^2.*exp(-theta*s)./H.^m;
den=tao*s.^3+(k2*h+1)*s.^2+(k2+k1*h)*s+k1;
T0=num./den;

T1 = T0.*T.^(n-1);

Mag = abs(T0);
Mag1 = abs(T1);
Max = max(Mag);
Max1 = max(Mag1(7:200));
if Max>1
    SS=1;
else
    SS=0.2*Max1;
end
y=p1*p2*p3*p4*SS;
end

