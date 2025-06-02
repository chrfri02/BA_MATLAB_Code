function f_mode = f_berechnet(N, l, S, C);

% Funktion zur Berechnung der Modenfrequenzen einer Metamaterialspule ueber die Birdcage-Regressionsformel
% dazu Verwendung mehrerer Regressionsmodelle aus den vorher aufgenommen Messwerten

% Angabe der Daten:
% Streifenanzahl N,
% Streifenlaenge l in mm, 
% Zellenabstand S in mm,
% Kapazitaet C in F

% Standardwerte, von denen ausgehend durch Variation eines Parameters die Messwerte aufgenommen wurden
N_Standard = 8;
l_Standard = 169;
S_Standard = 10;

% Koeffizienten, die durch Birdcage-Regression der Standardwerte ermittelt wurden 
L_Standard = 1.1731e-07;
alpha_Standard = 0.164600000000000;
k_Standard = 0.0837200000000000;

% Berechnung der zu den eingegebenen Parametern gehoerigen Koeffizienten
    % Induktivitaet L
        % absolute Aenderungen bezogen auf jeden Parameter
            % Streifenanzahl N
            a_L_N =  -4.0290e-08;
            b_L_N =   1.3133e-07;
            
            L_N = a_L_N*sqrt(1/N)+b_L_N;        % Modell sqrt(1/N)
            Delta_L_N = L_N-L_Standard;         % absolute Aenderung von L bezogen auf N abweichend vom Standardwert
    
            % Streifenlaenge l
            if l<=l_Standard
                a_L_l =   3.6327e-10;
                b_L_l =   5.5916e-08;
            else
                a_L_l =   7.4057e-10;
                b_L_l =  -7.8332e-09;
            end
    
            L_l = a_L_l*l+b_L_l;                % lineares Modell, aufgeteilt in zwei Intervalle
            Delta_L_l = L_l-L_Standard;         % absolute Aenderung von L bezogen auf l abweichend vom Standardwert
    
            % Zellenabstand S
            a_L_S =   4.4516e-09;
            b_L_S =   7.2842e-08;
    
            L_S = a_L_S*S+b_L_S;                % lineares Modell
            Delta_L_S = L_S-L_Standard;         % absolute Aenderung von L bezogen auf l abweichend vom Standardwert
        
        % Endwert fuer L als Summe des Standardwerts mit den absoluten Aenderungen bezogen auf jeden Spulenparameter
        L = L_Standard+Delta_L_N+Delta_L_l+Delta_L_S;   
    
    % Induktivitaetsverhaeltnis alpha
            % Streifenanzahl N
            a_alpha_N = 0.2246;
            b_alpha_N = 0.52;
            c_alpha_N = 0.1619;
        
            alpha_N = a_alpha_N*exp(-b_alpha_N*N)+c_alpha_N;    % exponentielles Modell mit negativem Vorzeichen im Exponenten
            Delta_alpha_N = alpha_N-alpha_Standard;             % absolute Aenderung von alpha bezogen auf N abweichend vom Standardwert
        
            % Streifenlaenge l
            a_alpha_l =    -0.00039;
            b_alpha_l =      0.2298;
        
            alpha_l = a_alpha_l*l+b_alpha_l;            % lineares Modell
            Delta_alpha_l = alpha_l-alpha_Standard;     % absolute Aenderung von alpha bezogen auf l abweichend vom Standardwert
        
            % Zellenabstand S
            if S <= S_Standard
               a_alpha_S =     0.00218;
               b_alpha_S =      0.1428;
            else
               a_alpha_S =     -0.0026;
               b_alpha_S =      0.1906;
            end
        
            alpha_S = a_alpha_S*S+b_alpha_S;            % lineares Modell, aufgeteilt in zwei Intervalle
            Delta_alpha_S = alpha_S-alpha_Standard;     % absolute Aenderung von alpha bezogen auf S abweichend vom Standardwert
    
        % Endwert fuer alpha als Summe des Standardwerts mit den absoluten Aenderungen bezogen auf jeden Spulenparameter
        alpha = alpha_Standard+Delta_alpha_N+Delta_alpha_l+Delta_alpha_S;
    
    % Koppelkonstante k
            % Streifenanzahl N
            a_k_N =     -0.1516;
            b_k_N =      0.1361;
    
            k_N = a_k_N*sqrt(1/N)+b_k_N;                % Modell sqrt(1/N)
            Delta_k_N = k_N-k_Standard;                 % absolute Aenderung von k bezogen auf N abweichend vom Standardwert
    
            % Streifenlaenge l
            if l <= l_Standard
                a_k_l =   -0.001509;
                b_k_l =      0.3388;
            else
                a_k_l =   5.567e-05;
                b_k_l =     0.07431;
            end
    
            k_l = a_k_l*l+b_k_l;                        % lineares Modell, aufgeteilt in zwei Intervalle
            Delta_k_l = k_l-k_Standard;                 % absolute Aenderung von k bezogen auf l abweichend vom Standardwert
    
            % Zellenabstand S
            if S <= S_Standard
                a_k_S =   -0.004156;
                b_k_S =      0.1253;
            else
                a_k_S =     0.00035;
                b_k_S =     0.08022;
            end
    
            k_S = a_k_S*S+b_k_S;                        % lineares Modell, aufgeteilt in zwei Intervalle
            Delta_k_S = k_S-k_Standard;                  % absolute Aenderung von k bezogen auf S abweichend vom Standardwert
    
        % Endwert fuer alpha als Summe des Standardwerts mit den absoluten Aenderungen bezogen auf jeden Spulenparameter
        k = k_Standard+Delta_k_N+Delta_k_l+Delta_k_S;

% Berechnung der Resonanzfrequenzen nach dem Birdcage-Modell

Omega_con = linspace(1,N-1);      % Modennummer Omega, kontinuierlich fuer den Graphen
Omega_dis = 1:N-1;                % Modennummer Omega, diskret fuer die Moden-Frequenzen

graph = 1./(2.*pi).*sqrt((1./(L.*C./2)).*(2.*sin(pi.*Omega_con./(2.*N)).*sin(pi.*Omega_con./(2.*N))./(alpha+2.*sin(pi.*Omega_con./(2.*N)).*sin(pi.*Omega_con./(2.*N)).*(1+2.*k.*cos(pi.*Omega_con./N)))));
f_mode = 1./(2.*pi).*sqrt((1./(L.*C./2)).*(2.*sin(pi.*Omega_dis./(2.*N)).*sin(pi.*Omega_dis./(2.*N))./(alpha+2.*sin(pi.*Omega_dis./(2.*N)).*sin(pi.*Omega_dis./(2.*N)).*(1+2.*k.*cos(pi.*Omega_dis./N)))));

% Ergebnisse plotten
fig1 = figure;
plot(Omega_con, graph, '-');
hold on;
grid on;
plot (Omega_dis, f_mode, 'x');
ylabel('Frequenz $$f$$ in Hz $$\rightarrow$$','Interpreter','latex','FontSize',18)
xlabel('Mode $$\Omega$$ $$\rightarrow$$','Interpreter','latex','FontSize',18)
legend('berechnete Frequenzen','Birdcage-Regressionkurve');
title(sprintf("berechnete Modenfrequenzen bei $$N=%d$$, $$l=%d$$ mm, $$S=%d$$ mm und $$C=%d$$ pF", N, l, S, C*10^12),'Interpreter','latex');
end