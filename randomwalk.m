clear
close all
delta_t=30;

part_conc=0.00001;

%vel_x=0.001;
vel_y=0.0000;
D=0.024;
L=3;
E_for=9*10^-7;
E_back=1*10^-8;
xmax=L; 
ymax=D/2;
x_conc=[0.6;1.2;2.4;3];


nbinx=30;
nbiny=20;

flow= [0 0;4 0.01;6 0.01; 30 0.25; 35 0.3;45 0.35; 90 0.36; 100 0.36];
conc_in= [0 0;4.99 0;5 250;6 250;6.01 0; 100 0] ;
% calcolo numero di particelle
dt=round(100*xmax/((max(flow(:,2))/1000)/(3.1415*(D^2)/4)*3/2)/nbinx/2)/100;
% dt=0.01;
tempo=0:dt:max(flow(:,1));
N=zeros(length(tempo),1);

for i=2:length(tempo)
    n_parc(i)=round(interp1(flow(:,1),flow(:,2),tempo(i))*...
        interp1(conc_in(:,1),conc_in(:,2),tempo(i))*dt/part_conc);
    N(i)=N(i-1)+n_parc(i);
end
x=zeros(max(N),1);
y=zeros(max(N),1);

NT=ceil(delta_t/dt);
conc=zeros(NT,length(x_conc)+1);


y=-ymax + 0.5*rand(max(N),1)*ymax;

Vx=linspace(0,xmax,nbinx+1);
Vy=linspace(-ymax,ymax,nbiny+1);


xsv(1)=x(1,1);
ysv(1)=y(1,1);

for n=1:NT
    t=n*(delta_t/NT);
    vel_x(n)=interp1(flow(:,1),flow(:,2),t)/1000/(3.1415*(D^2)/4);
    Re=vel_x(n)*D/(10^-6);
if Re<4200
    lam_check=1;
else
    lam_check=0;
end

    if lam_check==0
        for i = 1:N(n)
            r=randn(1,1);
            if r>0
                x(i,1)=x(i,1)+vel_x(n)*dt+r*(2*E_for*dt)^0.5;
            else
                x(i,1)=x(i,1)+vel_x(n)*dt+r*(2*E_back*dt)^0.5;
            end
        end
    else
        parfor i = 1:N(n)
            r=randn(1,1);
            if r>0
                x(i,1)=x(i,1)+3/2*vel_x(n)*(1.001-(y(i,1)/(D/2)).^2)*dt+r*(2*E_for*dt)^0.5;
            else
                x(i,1)=x(i,1)+3/2*vel_x(n)*(1.001-(y(i,1)/(D/2)).^2)*dt+r*(2*E_back*dt)^0.5;
            end
        end
    end
    
    parfor i = 1:N(n)
        r=randn(1,1);
        y(i,1)=y(i,1)+vel_y*dt+r*(10*(E_for+E_back)*dt)^0.5;
    end
    conc(n,1)=tempo(1,n);
    parfor i=1:N(n)
        if y(i,1)<-ymax
            y(i,1)=-2*ymax-y(i,1);
        else
            if y(i,1)>ymax
                y(i,1)=2*ymax-y(i,1);
            end
        end
    end
    for i=1:N(n)
        for k=1:length(x_conc)
            if (x(i,1)>=x_conc(k)) && (x(i,1)<(x_conc(k)+xmax/nbinx))
                conc(n,k+1)=conc(n,k+1)+part_conc/((xmax/nbinx)*3.1415*(D^2)/4)/1000;
            end
        end
    end
    
%     xsv(n+1)=x(1,1);
%     ysv(n+1)=y(1,1);
%     nx=histc(x,Vx);
%     ny=histc(y,Vy);
%     subplot(2,2,3);
%     bar(Vx,nx,'histc');
%     axis([0 xmax 0 max(N)/20]) 
%     xlabel('x'); ylabel('n');
%     title('particle conc. versus x')
%     subplot(2,2,2);
%     bar(Vy,ny,'histc');
%     xlabel('y'); ylabel('n');
%     title('particle conc. versus y')
%     axis([-ymax ymax 0 max(N)/4]) 
%     [count2d,x2d,y2d] = hist2d(x,y,nbinx,nbiny);
%     subplot(2,2,1)
%     plot(x,y,'ko',xsv,ysv,'r-o')
%     axis([0 xmax -ymax ymax])
%     xlabel('x'); ylabel('y'); title('particle position')
%     subplot(2,2,4)
%     image(x2d(1,:),y2d(:,1),count2d); 
%     colorbar;
%     axis([0 xmax -ymax ymax]);
%     set(gca,'YDir','normal');
%     title('2D concentration');
%     
%     pause(0.1)
   
end
figure
    for i=1:length(x_conc)
    plot(conc(:,1),conc(:,i+1))
    axis([0 delta_t 0 1.1*max(conc(:,2))]);
    xlabel('time'); ylabel('concentration');
    title('Concentration')
    hold on
    end
    writematrix(conc,'conc_lam_1.xls');
