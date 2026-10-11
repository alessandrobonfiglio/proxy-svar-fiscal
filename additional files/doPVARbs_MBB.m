function VARbs = doPVARbs_MBB(VAR, DAT, nboot, clevel, L)

% -------------------------------------------------------------
% Moving Block Bootstrap for Mertens & Ravn-style Proxy SVAR
% -------------------------------------------------------------

% Recompute reduced-form residuals (same as in doPVAR)
VARSLAGS = lagmatrix(VAR.vars,1:VAR.p);
VARSLAGS = VARSLAGS(VAR.p+1:end,:);
VARS     = VAR.vars(VAR.p+1:end,:);
res      = VARS - [VARSLAGS VAR.DET(VAR.p+1:end,:)]*VAR.bet;   % (VAR.T x VAR.n)

T = VAR.T;
N = VAR.n;
k = VAR.k;

TAXSHOCKS = VAR.TAXSHOCKS;   % aligned (VAR.T x VAR.k)

% -----------------------------
% Build moving blocks
% -----------------------------
J = ceil(T / L);

UBlocks = zeros(L, N, T-L+1);
MBlocks = zeros(L, k, T-L+1);

for ii = 1:T-L+1
    UBlocks(:,:,ii) = res(ii:ii+L-1,:);
    MBlocks(:,:,ii) = TAXSHOCKS(ii:ii+L-1,:);
end

% -----------------------------
% Position-in-block centering
% -----------------------------
centering  = zeros(L, N);
Mcentering = zeros(L, k);

for ii = 1:L
    centering(ii,:) = mean(res(ii:T-L+ii,:), 1);

    for jj = 1:k
        subM = TAXSHOCKS(ii:T-L+ii, jj);
        if any(subM ~= 0)
            Mcentering(ii,jj) = mean(subM(subM ~= 0));
        else
            Mcentering(ii,jj) = 0;
        end
    end
end

centering  = repmat(centering,  [J,1]); centering  = centering(1:T,:);
Mcentering = repmat(Mcentering, [J,1]); Mcentering = Mcentering(1:T,:);

% -----------------------------
% NEW: filters / diagnostics
% -----------------------------
% Your proxy is sparse -> require a minimum number of non-zero observations
minNZ = 20;      % try 10, 13, 15 (you had ~26 non-zero in the sample)
minAbsD11 = 1e-4; % try 1e-4 (looser) or 5e-3 (stricter)

discard_allzero = 0;
discard_weakNZ  = 0;
discard_failPVAR = 0;
discard_smallD11 = 0;

% -----------------------------
% Bootstrap loop
% -----------------------------
jj = 1;
while jj < nboot + 1

    % draw blocks
    idx = randi(T-L+1, J, 1);

    resb = zeros(J*L, N);
    taxb = zeros(J*L, k);

    for j = 1:J
        resb(1+L*(j-1):L*j,:) = UBlocks(:,:,idx(j));
        taxb(1+L*(j-1):L*j,:) = MBlocks(:,:,idx(j));
    end

    % IMPORTANT: trim to T
    resb = resb(1:T,:);
    taxb = taxb(1:T,:);

    % discard: all-zero proxy
    if all(taxb(:) == 0)
        discard_allzero = discard_allzero + 1;
        continue
    end

    % discard: too few non-zero proxy observations (weak instrument draw)
    if sum(taxb(:) ~= 0) < minNZ
        discard_weakNZ = discard_weakNZ + 1;
        continue
    end

    % center
    resb = resb - centering;
    for j = 1:k
        nz = (taxb(:,j) ~= 0);
        taxb(nz,j) = taxb(nz,j) - Mcentering(nz,j);
    end

    % --- FILTER 3: first-stage strength (rule-of-thumb: F > 10) ---
u_tax = resb(:,1);              % residual of tax revenue equation (adjust if tax is not col 1)
z     = taxb(:,1);              % proxy (k=1)

nz = (z ~= 0);
nEff = sum(nz);

if nEff < 10
    continue
end

X = [ones(nEff,1), z(nz)];
y = u_tax(nz);

b = X\y;
e = y - X*b;

SSR = sum(e.^2);
SST = sum((y - mean(y)).^2);

if SST <= 0
    continue
end

R2 = 1 - SSR/SST;
m  = 1;                         % number of instruments (k=1)
F  = (R2/(1-R2)) * ((nEff - m - 1)/m);

if ~isfinite(F) || F < 10
    continue
end


    % ---------------------------------------------------------
    % Generate bootstrap data using the SAME recursion as original
    % ---------------------------------------------------------
    varsb = zeros(VAR.p + T, N);
    varsb(1:VAR.p,:) = VAR.vars(1:VAR.p,:);

    for t = VAR.p+1 : VAR.p+T
        lvars = (varsb(t-1:-1:t-VAR.p,:))';
        varsb(t,:) = lvars(:)'*VAR.bet(1:VAR.p*N,:) ...
                   + VAR.DET(t,:)*VAR.bet(VAR.p*N+1:end,:) ...
                   + resb(t-VAR.p,:);
    end

    % Build bootstrap VAR object
    VARBS = VAR;
    VARBS.vars = varsb;
    VARBS.taxshocks = [VAR.TAXSHOCKS(1:VAR.p,:); taxb];  % pad first p

    % Re-estimate + identify (robust to failures)
    try
        VARBS = doPVAR(VARBS);
    catch
        discard_failPVAR = discard_failPVAR + 1;
        continue
    end

    % discard: unstable normalization for the tax shock (k=1 case)
    % Tax shock is j=1 when k=1 and tax variable is the first variable.
    if k == 1
        d11 = VARBS.D(1,1,1);
        if ~isfinite(d11) || abs(d11) < minAbsD11
            discard_smallD11 = discard_smallD11 + 1;
            continue
        end
    end

    % Store IRFs
    for j = 1:VAR.k
        irs = VARBS.irs(:,:,j);               % (irhor x N)
        VARbs.irs(:,jj,j) = irs(:);           % vectorize like original code
    end

    VARbs.irsg(:,jj)   = VARBS.irsg(:);
    VARbs.sigmaG(:,jj) = VARBS.sigmaG;
    VARbs.sigmaY(:,jj) = VARBS.sigmaY;
    VARbs.SigmaT(:,jj) = VARBS.SigmaT(:);
    VARbs.thetaY(:,jj) = VARBS.thetaY;
    VARbs.thetaG(:,jj) = VARBS.thetaG;
    VARbs.gammaT(:,jj) = VARBS.gammaT;
    VARbs.zetaT(:,jj)  = VARBS.zetaT;
    VARbs.zetaG(:,jj)  = VARBS.zetaG;
    VARbs.RM(:,jj)     = VARBS.RM(:);

    if VAR.k == 1
        VARbs.irsTRY(:,jj) = irs(:,1) - irs(:,3);
        VARbs.BS_taxshocks(:,jj) = VARBS.taxshocks;
        VARbs.BS_VARS(:,jj) = varsb(:);
    end

    jj = jj + 1;
end

% Save discard diagnostics (optional but helpful)
VARbs.discard_allzero  = discard_allzero;
VARbs.discard_weakNZ   = discard_weakNZ;
VARbs.discard_failPVAR = discard_failPVAR;
VARbs.discard_smallD11 = discard_smallD11;
VARbs.L = L;
VARbs.minNZ = minNZ;
VARbs.minAbsD11 = minAbsD11;

% -----------------------------
% Quantile bands (same style as original)
% -----------------------------
for j = 1:VAR.k
    VARbs.irsH(:,:,j) = reshape(quantile(VARbs.irs(:,:,j)', (1-clevel/100)/2), VAR.irhor, VAR.n);
    VARbs.irsL(:,:,j) = reshape(quantile(VARbs.irs(:,:,j)', 1-(1-clevel/100)/2), VAR.irhor, VAR.n);
end

if VAR.k == 1
    VARbs.irsTRYH = reshape(quantile(VARbs.irsTRY', (1-clevel/100)/2), VAR.irhor, 1);
    VARbs.irsTRYL = reshape(quantile(VARbs.irsTRY', 1-(1-clevel/100)/2), VAR.irhor, 1);
end

VARbs.irsgH = reshape(quantile(VARbs.irsg', (1-clevel/100)/2), VAR.irhor, VAR.n);
VARbs.irsgL = reshape(quantile(VARbs.irsg', 1-(1-clevel/100)/2), VAR.irhor, VAR.n);

VARbs.sigmaGci = [quantile(VARbs.sigmaG',(1-clevel/100)/2) quantile(VARbs.sigmaG',1-(1-clevel/100)/2)];
VARbs.sigmaYci = [quantile(VARbs.sigmaY',(1-clevel/100)/2) quantile(VARbs.sigmaY',1-(1-clevel/100)/2)];
VARbs.SigmaTci = [quantile(VARbs.SigmaT',(1-clevel/100)/2) quantile(VARbs.SigmaT',1-(1-clevel/100)/2)];
VARbs.thetaYci = [quantile(VARbs.thetaY',(1-clevel/100)/2)' quantile(VARbs.thetaY',1-(1-clevel/100)/2)'];
VARbs.thetaGci = [quantile(VARbs.thetaG',(1-clevel/100)/2)' quantile(VARbs.thetaG',1-(1-clevel/100)/2)'];
VARbs.gammaTci = [quantile(VARbs.gammaT',(1-clevel/100)/2)' quantile(VARbs.gammaT',1-(1-clevel/100)/2)'];
VARbs.zetaGci  = [quantile(VARbs.zetaG',(1-clevel/100)/2) quantile(VARbs.zetaG',1-(1-clevel/100)/2)];
VARbs.zetaTci  = [quantile(VARbs.zetaT',(1-clevel/100)/2)' quantile(VARbs.zetaT',1-(1-clevel/100)/2)'];
VARbs.RMci     = [quantile(VARbs.RM',(1-clevel/100)/2)' quantile(VARbs.RM',1-(1-clevel/100)/2)'];

end
