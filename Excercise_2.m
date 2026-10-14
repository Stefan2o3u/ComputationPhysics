clear;
N=17;          % beads
Nb=N-1;       % bonds
zeta=1;           % friction (unit)
kT=1;           % thermal energy (unit)
K=3;           % dimensionless spring constant
dt=0.01;     
M=50000;       % number of steps
randomForce=false %set true for assignment 3

%variables for the animation
nframes=50;
frames=floor(M/nframes);
rSaves= zeros(nframes, N, 3);    
tSaves=zeros(nframes, 1);
traj_end=zeros(M,1);                         

iframe = 0;
%Intitial condition
r= cumsum(randn(N,3), 1); %random starting position

rcm0=mean(r,1);
for step = 1:M
    

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
    if randomForce
        r=r+2*kT*zeta*dt*randn(N,3);
    end
    traj_end(step) = sum((r(N,:) - r(1,:)).^2);

    % save a frame
    if mod(step, nframes) == 0
        iframe= iframe + 1;
        rSaves(iframe,:,:) = r;
        tSaves(iframe)= step*dt;
    end
end
L = 0.6*Nb;                                  
for k= 1:iframe
    rk=squeeze(rSaves(k,:,:));              
    rk=rk-mean(rk,1);        %put it in frame of centre of mass            
    plot3(rk(:,1), rk(:,2), rk(:,3), 'o-', 'MarkerFaceColor','b');
    axis equal; grid on;
    xlim([-3 3]); ylim([-3 3]); zlim([-3 3]);
    xlabel('x'); ylabel('y'); zlabel('z');
    title(sprintf('t = %.2f', tSaves(k)));
    drawnow limitrate;
end


