%% Brownian dynamics of a Rouse chain - Exercise 4 (a and b)
clear; clc; close all;

%% Parameters (units: zeta = kT = l0 = 1)
N      = 17;          % beads
Nb     = N - 1;       % bonds
zeta   = 1;           % friction coefficient
kT     = 1;           % thermal energy
K      = 3;           % spring constant
dt     = 0.001;       % time step
M      = 2e7;         % number of steps (total time = M*dt = 20000)
nSave  = 10;          % store the bonds every nSave steps
nSamp  = M/nSave;     % number of stored samples
dtSamp = nSave*dt;    % time between stored samples

%% Transformation matrix B (q = B*u)
[j,p] = meshgrid(1:Nb, 1:Nb);
B = sqrt(2/(Nb+1)) * sin(j.*p*pi/(Nb+1));

%% Simulation
r = zeros(N,3);
%r(:,1)= 3*kT/K;
r(2:N,1) = cumsum(sqrt(kT/K)*randn(Nb,1));
noise = sqrt(2*kT*dt/zeta);                 % amplitude of the random displacement
U = zeros(nSamp, Nb);                       % x-projection of the bonds, u_j(t)

%optimized model
for step = 1:M
    bond = r(2:N,:) - r(1:N-1,:); % bond vectors (Nb x 3)

    F = zeros(N,3);                         % spring forces
    F(1:N-1,:) = F(1:N-1,:) + K*bond;       % bond i pulls bead i towards bead i+1
    F(2:N,:)   = F(2:N,:)   - K*bond;       % and bead i+1 towards bead i

    r = r + (dt/zeta)*F + noise*randn(N,3); % step foward + random force

    if mod(step, nSave) == 0                % save the bonds (X projection)
        U(step/nSave,:) = (r(2:N,1) - r(1:N-1,1))';
    end
end

%% Normal modes (part a) 
Q = U * B;                                  % q_p(t), one row per stored time
q2 = sum(Q.^2, 2);% |q|^2 at every step
u2 = sum(U.^2, 2);% |u|^2 at every step
fprintf('a) max | |q|^2 - |u|^2 | = %.2e\n', max(abs(q2 - u2)));
fprintf('   <q_p^2> should be kT/K = %.3f, found between %.3f and %.3f\n', ...
        kT/K, min(mean(Q.^2)), max(mean(Q.^2)));

%% Autocorrelation function of every mode, time-averaged 
maxLag  = 2000; % longest lag in samples (= 20 time units)
lagTime = (0:maxLag)' * dtSamp;
C = zeros(maxLag+1, Nb);
for s = 0:maxLag
    C(s+1,:) = mean(Q(1:end-s,:) .* Q(1+s:end,:), 1);
end
C = C ./ C(1,:);                            % divide by <q_p^2>, so C_p(0) = 1

%% Theory value and fit of ln C = -t/tau 
k = (1:Nb)' * pi/(Nb+1);
tauTheory = zeta ./ (2*K*(1 - cos(k)));
tauSim = zeros(Nb,1);
for ip = 1:Nb
    below = find(C(:,ip) < 0.2, 1);         % fit only while C_p is above 0.2 to reduce noise
    if isempty(below), below = maxLag+1; end
    last = below - 1;
    fit = polyfit(lagTime(1:last), log(C(1:last,ip)), 1);% fit curve
    tauSim(ip) = -1/fit(1);
end

fprintf('\n  p   tau_theory   tau_sim   deviation\n');
for ip = 1:Nb
    fprintf('%3d %10.4f %10.4f %8.1f%%\n', ip, tauTheory(ip), tauSim(ip), ...
            100*(tauSim(ip)/tauTheory(ip) - 1));
end

%% Plots
figure;
modes = [1 5 8 12 16];
for i = 1:numel(modes)
    ip = modes(i);
    subplot(2,3,i);
    semilogy(lagTime, C(:,ip), 'b.', lagTime, exp(-lagTime/tauTheory(ip)), 'r-', 'LineWidth', 1.5);
    xlim([0 min(4*tauTheory(ip), lagTime(end))]);
    ylim([0.02 1]);
    xlabel('t'); ylabel('C_p(t)'); title(sprintf('p = %d', ip));
    if i == 1, legend('simulation', 'e^{-t/\tau_p} (theory)'); end
end
subplot(2,3,6);
semilogy(1:Nb, tauTheory, 'r-', 1:Nb, tauSim, 'bo', 'LineWidth', 1.5);
xlabel('p'); ylabel('\tau_p'); legend('theory', 'simulation');