function GVS_model_IMM_clean()

    function Phi=Strain_Modes(s)
        Phi=[[0;0;0;1;0;0],...
             [0;0;0;0;1;0],...
             [0;0;0;0;0;1],...
             [0;0;0;s;0;0],...
             [0;0;0;0;s;0],...
             [0;0;0;0;0;s],...
             [0;0;0;s^2;0;0],...
             [0;0;0;0;s^2;0],...
            [0;0;0;0;0;s^2]];
    end

    function xihat=hat6(xi)
        xihat=[    0, -xi(6), xi(5), xi(1);
               xi(6),      0,-xi(4), xi(2);
              -xi(5),  xi(4),     0, xi(3);
              0,0,0,0];
    end

    function uhat=hat3(u)
        uhat=[    0, -u(3), u(2);
               u(3),      0,-u(1);
              -u(2),  u(1),     0];
    end

    function adxi=ad(xi)
        v=xi(1:3); vhat=hat3(v); u=xi(4:6); uhat=hat3(u);
        adxi=[ uhat, vhat; 
               zeros(3,3), uhat];
    end

    function adbarW=adbar(W)
        f=W(1:3); fhat=hat3(f); m=W(4:6); mhat=hat3(m);
        adbarW=-[ zeros(3,3), fhat; 
               fhat, mhat];
    end

    function x_prime=DynamicsODE_with_C(s,x)

        T=reshape(x(1:16),4,4); 
        x(1:16)=[];
        J=reshape(x(1:6*n),6,n); 
        x(1:6*n)=[];
        A=reshape(x(1:6*n),6,n);

        Phi=Strain_Modes(s);
        xi=Phi*q+xi_star;
        xidot=Phi*qd;
        adxi=ad(xi);
        adxidot=ad(xidot);
        V=J*qd;
        adV=ad(V);
        adbarMV=adbar(Md*V);

        R=T(1:3,1:3);
        wb=[R.'*a_grav*rho*Area;0;0;0];

        T_prime=T*hat6(xi);
        J_prime=Phi-adxi*J;
        M_prime=J.'*Md*J;
        A_prime=adV*Phi-adxidot*J-adxi*A;
        C_prime=1/2*J.'*(Md*A-adV.'*Md*J+adbarMV*J);
        Q_prime=J.'*wb;
        PE_prime=1/2*(xi-xi_star).'*Kd*(xi-xi_star)-T(1:3,4).'*a_grav*rho*Area;
        KE_prime=1/2*V.'*Md*V;

        x_prime=[   reshape(T_prime,16,1);
                    reshape(J_prime,6*n,1);
                    reshape(A_prime,6*n,1);
                    reshape(M_prime,n*n,1);
                    reshape(C_prime,n*n,1);
                    reshape(Q_prime,n,1);
                    PE_prime;KE_prime];
    end

    function x_prime=DynamicsODE(s,x)

        T=reshape(x(1:16),4,4); 
        x(1:16)=[];
        J=reshape(x(1:6*n),6,n); 
        x(1:6*n)=[];
        a=x(1:6);

        Phi=Strain_Modes(s);
        xi=Phi*q+xi_star;
        xidot=Phi*qd;
        adxi=ad(xi);
        V=J*qd;
        adV=ad(V);

        R=T(1:3,1:3);
        wb=[R.'*a_grav*rho*Area;0;0;0];

        T_prime=T*hat6(xi);
        J_prime=Phi-adxi*J;
        M_prime=J.'*Md*J;
        a_prime=adV*xidot-adxi*a;
        c_prime=J.'*(Md*a-adV.'*Md*V);
        Q_prime=J.'*wb;
        PE_prime=1/2*(xi-xi_star).'*Kd*(xi-xi_star)-T(1:3,4).'*a_grav*rho*Area;
        KE_prime=1/2*V.'*Md*V;

        x_prime=[   reshape(T_prime,16,1);
                    reshape(J_prime,6*n,1);
                    a_prime;
                    reshape(M_prime,n*n,1);
                    c_prime;Q_prime;PE_prime;KE_prime];
    end

    function [t,Y]=RK4(ODE,t1,t2,y1,N)
        h1=(t2-t1)/N;
        Y=zeros(length(y1),n);
        t=linspace(t1,t2,N+1);
        Y(:,1)=y1;
        for i=1:N
            k1=ODE(t(i),Y(:,i));
            k2=ODE(t(i)+h1/2,Y(:,i)+k1*h1/2);
            k3=ODE(t(i)+h1/2,Y(:,i)+k2*h1/2);
            k4=ODE(t(i)+h1,Y(:,i)+k3*h1);
            Y(:,i+1)=Y(:,i)+h1/6*(k1+2*k2+2*k3+k4);
        end
        Y=Y.';
    end

    function [qkp1,qdkp1,tkp1,X,j]=IMM_step(qk,qdk,tk,h)
        
        q=qk+h/2*qdk;
        qd=qdk;
        num_iter=20;
        for j=1:num_iter
            %[S,X]=ode45(@DynamicsODE,0:.02:L,x0);
            [S,X]=RK4(@DynamicsODE,0,L,x0,30);
            x=X(end,:);
            Tend=reshape(x(1:16),4,4)  ;
            x(1:16)=[];
            J=reshape(x(1:6*n),6,n) ;   
            x(1:6*n)=[];
            a=x(1:6).';                 
            x(1:6)=[];
            M=reshape(x(1:n*n),n,n);
            x(1:n*n)=[];
            c=x(1:n).';                 
            x(1:n)=[];
            Q=x(1:n).';              


            H=2/h*M+h/2*K+B;
            b=2/h*M*qdk-c-K*qk+Q;
            
            qd_new=H\b;
            q_new=qk+h/2*qd_new;
            if norm(q_new-q)/n<.0001 &&  norm(qd_new-qd)/n*h<.0001
                q=q_new;
                qd=qd_new;
                break
            end
                q=q_new;
                qd=qd_new;
        end
        
        qkp1=2*q-qk;
        qdkp1=2*qd-qdk;
        tkp1=tk+h;
    end


    function [qkp1,qdkp1,tkp1,X,j]=IMM_step_C(qk,qdk,tk,h)
        
        q=qk+h/2*qdk;
        qd=qdk;
        num_iter=20;
        for j=1:num_iter
            %[S,X]=ode45(@DynamicsODE_with_C,0:.02:L,x0);
            [S,X]=RK4(@DynamicsODE_with_C,0,L,x0,30);
            x=X(end,:);
            Tend=reshape(x(1:16),4,4);  x(1:16)=[];
            J=reshape(x(1:6*n),6,n);    x(1:6*n)=[];
            A=reshape(x(1:6*n),6,n);    x(1:6*n)=[];
            M=reshape(x(1:n*n),n,n);    x(1:n*n)=[];
            C=reshape(x(1:n*n),n,n);    x(1:n*n)=[];
            Q=x(1:n).';                   %x(1:n)=[];

            H=2/h*M+2*C+h/2*K+B;
            b=2/h*M*qdk+C*qd-K*qk+Q;

            %H=2/h*M+h/2*K+B;
            %b=2/h*M*qdk-C*qd-K*qk+Q;
            
            qd_new=H\b;
            q_new=qk+h/2*qd_new;
            if norm(q_new-q)/n<.0001 &&  norm(qd_new-qd)/n*h<.0001
                q=q_new;
                qd=qd_new;
                break
            end
                q=q_new;
                qd=qd_new;
        end
        
        qkp1=2*q-qk;
        qdkp1=2*qd-qdk;
        tkp1=tk+h;
    end

clear all
close all
r=1e-3; %meters
L=1; % meters
Area=pi*r^2;
Ixx=pi*r^4/4;
Iyy=pi*r^4/4;
Izz=pi*r^4/4;
E=60e9; %Pa
nu=0.3;
G=E/(2*(1+nu)); %Pa
rho=6500; %kg/m^3
xi_star=[0;0;1;1;0;0];

global K B Md q qd f_grav
n=size(Strain_Modes(0),2);
Kd=diag([Area*G,Area*G,Area*E,E*Ixx,E*Iyy,G*Izz]);
alpha=0;
Md=diag([Area*rho,Area*rho,Area*rho,rho*Ixx,rho*Iyy,rho*Izz]);
K_integrand=@(s,y) reshape(Strain_Modes(s).'*Kd*Strain_Modes(s),n*n,1);
[~,Karray]=ode45(K_integrand,[0 L],zeros(n,n));
K=reshape(Karray(end,:),n,n)
B=alpha*K;
a_grav=[-9.81;0;0];
q=zeros(n,1);
qd=zeros(n,1);

global x0

% Use these for the version with C
% T0=eye(4);
% J0=zeros(6,n);
% A0=zeros(6,n);
% M0=zeros(n,n);
% C0=zeros(n,n);
% Q0=zeros(n,1);
% 
% x0=[reshape(T0,16,1);
%     reshape(J0,6*n,1);
%     reshape(A0,6*n,1);
%     reshape(M0,n*n,1);
%     reshape(C0,n*n,1);
%     reshape(Q0,n,1);0;0];


%Use these for the version with c
T0=eye(4);
J0=zeros(6,n);
a0=zeros(6,1);
M0=zeros(n,n);
c0=zeros(n,1);
Q0=zeros(n,1);
x0=[reshape(T0,16,1);
    reshape(J0,6*n,1);
    a0;
    reshape(M0,n*n,1);
    c0;Q0;0;0];


q0=zeros(n,1);
qd0=zeros(n,1);
y0=[q0;qd0];

qk=q0;
qdk=qd0;
tk=0; 
h=.02;

N_steps=200;
q_list=zeros(n,N_steps+1);
qd_list=zeros(n,N_steps+1);
t_list=linspace(0,h*N_steps,N_steps+1);

q_list(:,1)=q0;
qd_list(:,1)=qd0;

totalj=0;
for k=1:N_steps
        
        [qkp1,qdkp1,tkp1,X,j]=IMM_step(qk,qdk,tk,h);
        
        q_list(:,k+1)=qkp1;
        qd_list(:,k+1)=qdkp1;
        PE_list(k)=X(end,end-1);
        KE_list(k)=X(end,end);

        p=X(:,13:15);
        plot3(p(:,1),p(:,2),p(:,3));
        grid on
        view([1 1 .5])
        xlim([-1 1])
        ylim([-1 1])
        zlim([0 1])
        hold off
        daspect([1 1 1])  
        drawnow
        qk=qkp1;
        qdk=qdkp1;
        tk=tkp1;

        totalj=totalj+j;
end

KE_list

figure
plot(t_list,q_list)
xlabel('Time')
ylabel('q')

figure
plot(t_list(2:end),KE_list,t_list(2:end),PE_list,t_list(2:end),KE_list+PE_list)
xlabel('Time')
ylabel('Energy')
legend('Kinetic','Potential','Total')



end