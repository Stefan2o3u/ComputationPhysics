N=17;          % beads
Nb=N-1;       % bonds
zeta=1;           % friction (unit)
kT=1;           % thermal energy (unit)
K=3;           % dimensionless spring constant
dt=0.001;     
M=50000;       % number of steps
randomForce=true %set true for assignment 3
meanBond=zeros(M,1);

%variables for the animation
nframes=50;
frames=floor(M/nframes);
rSaves= zeros(nframes, N, 3);    
tSaves=zeros(nframes, 1);
bondTimeStep=zeros(M,Nb);

% EX4:Transformation matrix B (Nb x Nb)
[j,p] = meshgrid(1:Nb, 1:Nb);
B = sqrt(2/(Nb+1)) * sin(j.*p*pi/(Nb+1));
fprintf('max |B*B - I| = %.2e\n', max(max(abs(B*B - eye(Nb)))));

iframe = 0;
%Intitial condition
%r= cumsum(randn(N,3), 1); %random starting position
%r = [(0:N-1)', zeros(N,3-1)];
r = zeros(N, 3);
r(:,1)= 3*kT/K
rcm0=mean(r,1);
U = zeros(M,Nb);
for step = 1:M
    b=r(2:N,:)-r(1:N-1,:);

    %SpringForce
    F = zeros(N,3);

    % bead 1
    F(1,:) = K*( r(2,:) - r(1,:) );
    % inner beads 
    for i = 2:N-1
        F(i,:) = K*( r(i+1,:) - 2*r(i,:) + r(i-1,:) );
    end
    % bead N
    F(N,:) = -K*( r(N,:) - r(N-1,:) );

    %Step forward
    r=r+(dt/zeta)*F;
    %Assignment 3
    r=r+sqrt(2*kT*dt/zeta)*randn(N,3);

    %Assignment 4
    U(step,:) = (r(2:N,1) - r(1:N-1,1))';   % bonds, X projection

    bondTimeStep(step,:)=(sum(b.^2, 2));
    meanBond(step)=mean(sum(b.^2, 2));

    % save a frame
    if mod(step, nframes)==1
        iframe= iframe + 1;
        rSaves(iframe,:,:) = r;
        tSaves(iframe)= step*dt;
    end
end
t = (1:M)*dt;
figure; plot(t, meanBond); xlabel('t'); ylabel('<l^2>_{bonds}')
figure; plot(t, bondTimeStep(:,1)); xlabel('t'); ylabel('l_1^2{bonds}')
timeAveraged=3*sqrt(pi)/8*(2/K)^(5/2)*4*pi*(K/(2*kT*pi))^(3/2)

L = 0.2*Nb;                                  
for k= 1:iframe
    rk=squeeze(rSaves(k,:,:));              
    rk=rk-mean(rk,1);        %put it in frame of centre of mass            
    plot3(rk(:,1), rk(:,2), rk(:,3), 'o-', 'MarkerFaceColor','b');
    axis equal; grid on;
    xlim([-L L]); ylim([-L L]); zlim([-L L]);
    xlabel('x'); ylabel('y'); zlabel('z');
    title(sprintf('t = %.2f', tSaves(k)));
    drawnow limitrate;
end

% Exercise 4
Q = U * B;
q2 = sum(Q.^2, 2);% |q|^2 at every step
u2 = sum(U.^2, 2);% |u|^2 at every step

fprintf('max | |q|^2 - |u|^2 | = %.2e\n', max(abs(q2-u2)));
fprintf('max relative difference = %.2e\n', max(abs(q2-u2)./u2));

t = (1:M)*dt;
figure;
subplot(2,1,1);
plot(t,u2,'b-',t,q2,'r--'); legend('|u|^2','|q|^2');
xlabel('t'); ylabel('squared length');
subplot(2,1,2);
semilogy(t, abs(q2-u2)); xlabel('t'); ylabel('| |q|^2 - |u|^2 |');