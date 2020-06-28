%Flexural isostasy
%  Te = Elastic thickness
% ero = load
% dx = x spacing
% dy = y spacing
% pc = load density
% buffer = edge buffer size (in grid node units);
function [w] = IsoFlex(ero,Te,dx,dy,pc,buffer)

[m,n] = size(ero);


continuous = true; %assume load is continuous over buffer zone or make load zero over buffer zone?
E = 100e9;
g = 9.81;
v = .25;
pm = 3300;
%buffer = 0;%100;
bufferx = 0;
buffery = 0;


pw = 00;%%% water density
if mod(n,2) == 1
    bufferx = bufferx+1;
end
if mod(m,2)==1
    buffery = buffery+1;
end
D = E*(Te)^3/(12*(1-v^2));

if continuous
    appndl = zeros(m,buffer+bufferx);
    for i = 1:buffer +bufferx
       appndl(:,i) =  ero(:,1);
    end
    appndr = zeros(m,buffer);
    for i = 1:buffer 
       appndr(:,i) =  ero(:,end);
    end
    ero = [appndl ero appndr];
    appndu = zeros(buffer+buffery,n+2*buffer +bufferx);
    for i = 1:buffer+buffery
        appndu(i,:) =  ero(1,:);
    end
    appndd = zeros(buffer,n+2*buffer+bufferx);
   for i = 1:buffer
        appndd(i,:) =  ero(end,:);
    end
    ero = [appndu ;ero ;appndd];
else
    appndlr = zeros(m,(buffer+bufferx));
    appndud = zeros((buffer+buffery),2*(buffer)+n+bufferx);
    ero = [appndlr ero appndlr(:,1:end-bufferx)];
    ero = [appndud;ero;appndud(1:end-buffery,:)];
end



[m,n] = size(ero);
Ly = m*dy;
Lx = n*dx;
h = fft2(ero);
k = zeros(m/2+1,n/2+1);
k(1,1) = pc/(pm-pw);%(1+(D/(g*(pm-pc))*(2*pi*(1)/(Ly)).^4));
for i = 2:ceil((m-1)/2)+1
      k(i,1) =  ((pc)/(pm-pw))*1/(1+(D/(g*(pm-pw))*(2*pi*(i-1)/(Ly)).^4));
end
for j = 2:ceil((n-1)/2)+1
      k(1,j) = (((pc)/(pm-pw))*1/(1+(D/(g*(pm-pw))*(2*pi*(j-1)/(Lx)).^4)));
end
for i = 2:ceil((m-1)/2)+1
    for j = 2:ceil((n-1)/2)+1
        k(i,j) = (((pc)/(pm-pw)).*1/(1+(D/(g*(pm-pw))*(2*pi*sqrt(((i-1)/(Ly))^2+((j-1)/(Lx))^2)).^4)));
    end
end
k = [k (fliplr(k(:,2:end-1)))];
k = [k;(flipud(k(2:end-1,:)))];
w = k.*h;
w = ifft2(w);
w = w(buffer+1:end-buffer-bufferx,buffer+1:end-buffer-buffery);
end

