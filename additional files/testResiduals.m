function residualTests = testResiduals(VAR, maxLag)
% Test di diagnostica sui residui del modello VAR
% Input:
%   VAR - struttura contenente i risultati del VAR
%   maxLag - numero massimo di ritardi da testare (default: 12)
%
% Output:
%   residualTests - struttura con i risultati dei test

if nargin < 2
    maxLag = 12; %se non viene specificato un numero di lag viene automaticamente impostato 12
end

% Calcola i residui
VARSLAGS = lagmatrix(VAR.vars, 1:VAR.p);
VARSLAGS = VARSLAGS(VAR.p+1:end, :);
VARS = VAR.vars(VAR.p+1:end, :);
res = VARS - [VARSLAGS VAR.DET(VAR.p+1:end, :)] * VAR.bet;
%sottrae manualmenete dai dati originali il valore predetto dal modello
%usando le matrici dei ritardi (VARLAGS), le costanti/trend (DET) e i
%coefficienti stimati (BET)
[T, n] = size(res); %estrae le dimensioni

% Salva i residui nella struttura output
residualTests.residuals = res;
residualTests.T = T;
residualTests.n = n;

%% 1. TEST DI LJUNG-BOX (multivariato)
fprintf('\n=== TEST DI LJUNG-BOX (Multivariato) ===\n');
% Inizializza i vettori per statistica Q, p-value e i gradi di libertà
residualTests.LjungBox.Q = zeros(maxLag, 1); 
residualTests.LjungBox.pvalue = zeros(maxLag, 1);
residualTests.LjungBox.df = zeros(maxLag, 1);

for lag = 1:maxLag
  
    Q = 0;
    for h = 1:lag
       
        C0 = (res' * res) / T;
        Ch = (res(1:end-h, :)' * res(h+1:end, :)) / T;
        Q = Q + trace(Ch' / C0 * Ch / C0) / (T - h);
    end
    Q = T * (T + 2) * Q;
    
    % Se il numero di lag è inferiore all'ordine del VAR, il test non è definibile.
    df = n^2 * (lag - VAR.p);
    % se i gradi di libertà sono positivi, calcola il p-value usando la
    % distribuzione chi-quadro
    if df > 0
        pvalue = 1 - chi2cdf(Q, df);
        residualTests.LjungBox.Q(lag) = Q;
        residualTests.LjungBox.pvalue(lag) = pvalue;
        residualTests.LjungBox.df(lag) = df;
        
        if lag <= 4 || lag == maxLag
            fprintf('Lag %2d: Q = %8.4f, df = %3d, p-value = %6.4f %s\n', ...
                lag, Q, df, pvalue, getSignificance(pvalue));
        end
    end
end

%% 2. TEST DI BREUSCH-GODFREY (per ogni equazione)
% Viene svolto equazione per equazione per cercare autocorrelazione specifica
fprintf('\n=== TEST DI BREUSCH-GODFREY ===\n');
% Inizializza i contenitori dei risultati per ogni equazione
residualTests.BreuschGodfrey.LM = zeros(n, 1);
residualTests.BreuschGodfrey.pvalue = zeros(n, 1);
residualTests.BreuschGodfrey.df = zeros(n, 1);

% Inizia il ciclo sulle n variabili del sistema
for eq = 1:n
    % Test LM per l'equazione eq
    lagBG = min(4, floor(T/5)); % Numero di lag per il test
    
    % Regressione dei residui sui ritardi
    resLags = lagmatrix(res(:, eq), 1:lagBG);
    resLags = resLags(lagBG+1:end, :);
    resEq = res(lagBG+1:end, eq); %allinea i residui attuali dell'equazione allineandoli alla matrice dei ritardi appena creata
    XorigLags = [VARSLAGS(lagBG+1:end, :) VAR.DET(VAR.p+lagBG+1:end, :)]; % recupera i dati originali usati nel VAR (variabili ritardate e costanti)
    
    % Regressione ausiliaria
    X_aux = [XorigLags resLags];
    beta_aux = X_aux \ resEq; % cerca di spiegare i residui attuali 'resEq' usando sia i vecchi dati che i residui passati
    res_aux = resEq - X_aux * beta_aux; % calcola l'errore rimanente. se i ritardi dei residui sono utili a spiegare 'resEq', questo nuovo errore sarà molto più piccolo dell'originale
    
    % Statistica LM
    SSR_res = sum(res_aux.^2); %calcola somma dei quadrati degli errori della nuova regressione
    SSR_orig = sum(resEq.^2); % somma dei quadrati degli errori originali
    T_eff = length(resEq); 
    LM = T_eff * (1 - SSR_res / SSR_orig); % formula del test
    
    df_BG = lagBG;
    pvalue_BG = 1 - chi2cdf(LM, df_BG); % p-value
    
    residualTests.BreuschGodfrey.LM(eq) = LM;
    residualTests.BreuschGodfrey.pvalue(eq) = pvalue_BG;
    residualTests.BreuschGodfrey.df(eq) = df_BG;
    
    fprintf('Equazione %d: LM = %8.4f, df = %2d, p-value = %6.4f %s\n', ...
        eq, LM, df_BG, pvalue_BG, getSignificance(pvalue_BG));
end

%% 3. TEST PORTMANTEAU
fprintf('\n=== TEST PORTMANTEAU ===\n');
residualTests.Portmanteau.Q = zeros(maxLag, 1);
residualTests.Portmanteau.pvalue = zeros(maxLag, 1);
residualTests.Portmanteau.df = zeros(maxLag, 1);
% il test inizia dal lag p+1. 
for lag = VAR.p+1:maxLag
    % Statistica di portmanteau
    Qport = 0;
    C0_inv = inv((res' * res) / T); 
    % Calcola l'inversa della matrice di var-cov dei residui al tempo zero.
    % normalizza i residui prima di misurarne la correlazione
    % inizia un ciclo interno per sommare i contributi di ogni ritardo fino
    % a quello corrente
    for h = 1:lag
        Ch = (res(1:end-h, :)' * res(h+1:end, :)) / T; % matrice di autocovarianza al ritardo h
        Qport = Qport + trace(Ch' * C0_inv * Ch * C0_inv) / (T - h); 
    end
    Qport = T^2 * Qport;
    
    df_port = n^2 * (lag - VAR.p); % calcola i gradi di libertà nel test. dipendono dal numero di variabili al quadrato e dalla distanza tra il lag attuale e l'ordine del var
    pvalue_port = 1 - chi2cdf(Qport, df_port); % calcola il p-value usando distribuzione chi-quadro
    
    residualTests.Portmanteau.Q(lag) = Qport;
    residualTests.Portmanteau.pvalue(lag) = pvalue_port;
    residualTests.Portmanteau.df(lag) = df_port;
    %qua stampa risultati solo per il primo e l'ultimo lag
    if lag == VAR.p+1 || lag == maxLag
        fprintf('Lag %2d: Q = %8.4f, df = %3d, p-value = %6.4f %s\n', ...
            lag, Qport, df_port, pvalue_port, getSignificance(pvalue_port));
    end
end

%% 4. AUTOCORRELAZIONI DEI RESIDUI
fprintf('\n=== AUTOCORRELAZIONI DEI RESIDUI ===\n');

% Determina il numero massimo di lag per ACF
maxLagACF = min(maxLag, floor(T/4));

% Inizializza le matrici con la dimensione corretta
residualTests.ACF.values = zeros(n, maxLagACF);
residualTests.ACF.se = zeros(n, maxLagACF);

for eq = 1:n
    for lag = 1:maxLagACF
        % Calcola autocorrelazione
        acf = corr(res(1:end-lag, eq), res(lag+1:end, eq));
        residualTests.ACF.values(eq, lag) = acf;
        
        % Deviazione standard approssimata
        se = 1/sqrt(T);
        residualTests.ACF.se(eq, lag) = se;
    end
    
    % Mostra solo i primi 4 lag
    numLagsToShow = min(4, maxLagACF);
    fprintf('Equazione %d - ACF(1:%d): [', eq, numLagsToShow);
    for lag = 1:numLagsToShow
        fprintf('%6.3f', residualTests.ACF.values(eq, lag));
        if lag < numLagsToShow
            fprintf(' ');
        end
    end
    fprintf(']\n');
end

% Salva il numero di lag calcolati
residualTests.ACF.maxLag = maxLagACF;

%% 5. PLOT DEI RESIDUI E ACF
figure('Name', 'Diagnostica Residui');

for eq = 1:n
    % Plot dei residui
    subplot(n, 3, (eq-1)*3 + 1);
    plot(res(:, eq), 'LineWidth', 1);
    hold on;
    yline(0, 'r--');
    yline(2*std(res(:, eq)), 'k:');
    yline(-2*std(res(:, eq)), 'k:');
    title(sprintf('Residui - Eq %d', eq));
    xlabel('Osservazioni');
    grid on;
    
    % ACF
    subplot(n, 3, (eq-1)*3 + 2);
    maxLagPlot = min(20, maxLagACF);
    acfVals = residualTests.ACF.values(eq, 1:maxLagPlot);
    bar(1:maxLagPlot, acfVals);
    hold on;
    yline(1.96/sqrt(T), 'r--');
    yline(-1.96/sqrt(T), 'r--');
    title(sprintf('ACF - Eq %d', eq));
    xlabel('Lag');
    ylabel('ACF');
    grid on;
    
    % Q-Q plot
    subplot(n, 3, (eq-1)*3 + 3);
    qqplot(res(:, eq));
    title(sprintf('Q-Q Plot - Eq %d', eq));
end

end

function sigStr = getSignificance(pval)
    if pval < 0.01
        sigStr = '***';
    elseif pval < 0.05
        sigStr = '**';
    elseif pval < 0.10
        sigStr = '*';
    else
        sigStr = '';
    end
end