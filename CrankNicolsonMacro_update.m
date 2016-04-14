
function [U,A,B] = CrankNicolsonMacro_update(intX,T,dt,D,l0)

% To solve the equation du/dt  = Dd^2u/dx^2 + f(u)
% discretize using Crank Nicolson and obtain the 
% system AU^{n+1} = BU^n + f^n, and so solving for 
% U^{n+1} gives us: U^{n+1} = A^{-1}(BU^n + \deltaT * f^n)
% To minimize computation, we invert A, a constant matrix, only once. 

% init
a = intX(1);
dx = intX(2)-a;
tao = D*(dt/2)/dx^2;

%% Create A & B
A = -tao*diag(ones(1,length(intX)-1),-1) - tao*diag(ones(1,length(intX)-1),1) ;
A(logical(eye(size(A)))) = 1+2*tao;
A(1,1) = 1+tao; A(1,2) = -tao;          % Neumann BC: u_0 = u_1, u_N+1 = u_N
A(end:end-1) = -tao; A(end,end) = 1+tao;

% Create B
B = tao*diag(ones(1,length(intX)-1),-1) + tao*diag(ones(1,length(intX)-1),1) ;
B(logical(eye(size(B)))) = 1-2*tao;
B(1,1) = 1-tao; B(1,2) = tao;
B(end:end-1) = tao; B(end,end) = 1-tao;

%% Solve for U^{n+1}
U = zeros(length(intX),length(T)); % solution matrix
%  U(:,1) = normpdf(intX,0,1); % initial condition
%  U(:,1) = U(:,1)/(sum(U(:,1)*dx));
nL0 = floor(l0/dx+.5);
nIndex0 =  floor( (0-a)/dx+.5) + 1;
U( (-nL0:nL0) + nIndex0, 1) = 1;
U(:,1) = U(:,1)/sum((U(:,1)*dx)); %normalize IC
f = @(y) 0*y.*(1-y);
A_inv = inv(A);


% figure;
for j = 2:length(T)
    % Strang splitting
    % ----------------
    % A) ∂_t u = f(u) over Δt/2
    U1 = U(:,j-1) + dt/2*f(U(:,j-1));
    % B) ∂_t u = ∂_xx u over Δt
    %U2 = A\(B*U1);
    U2 = A_inv*(B*U1);
    % C) ∂_t u = f(u) over Δt/2
    U(:,j) = U2 + dt/2*f(U2);
    % one shot
    % --------
    %U(:,j) = A\(B*U(:,j-1) + dt*f(U(:,j-1)));
%      plot(intX,U(:,j)); pause(.1); grid on;
    U(:,j) = U(:,j)/(sum(U(:,j))*dx); % normalize
end


% end







