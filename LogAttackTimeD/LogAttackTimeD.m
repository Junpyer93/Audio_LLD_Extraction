% function [LogAttackTime] = mp7DLogAttackTime(envelop_bp, threshold_percent)
%
% estimate the attack start and end time of a signal envelop
% and 
% compute the log-attack-time
%
% INPUTS
% ======
% - envelop_bp        : energy envelope (first column: time [second] | second column: value)
% - threshold_percent : percentage of maximum signal energy applied in
% order to determine start time (Solitamente per convenzione inserisco 2%)
% 
% OUTPUTS
% =======
% - LogAttackTime               : log-attack-time
%
% Target:   MP7-XM version
% Author:   CUIDADO/IRCAM/ G. Peeters 
% LastEdit: 2001/03/12
%

function [LogAttackTime, StAtt, StPAtt] = LogAttackTimeD(envelop_bp, threshold_percent)
  
  time_v   = envelop_bp(:,1); %tutti i valori della prima colonna
  energy_v = envelop_bp(:,2); %tutti i valori della seconda colonna, usa la funzione envelope(X)
  
  % === informative
  [stopattack.value, stopattack.pos] = max(energy_v); %trovo il punto in cui si ha il massimo valore d'energia(valore, posizione nel vettore)
  threshold = stopattack.value * threshold_percent/100; %Dato il treshold in percentuale, trovo il valore associato
  tmp       = find(energy_v > threshold); %trova l'indice lineare dei valori del vettore di energia che sono più grandi della treshold
  disp(tmp);
  startattack.pos = tmp(1);
  disp(tmp(1));
  
  if (startattack.pos == stopattack.pos)
      startattack.pos = startattack.pos - 1;
  end
  
  % === normative
  
  disp(threshold);
  disp(startattack.pos);
  disp(stopattack.pos);
  disp(stopattack.value);
  disp(time_v(stopattack.pos));
  %disp(time_v(startattack.pos));
  
  % Il seguente if è stato imposto affinchè si preveda il caso in cui
  % startattack.pos =0. In quel caso la dimensione di time_v è inesistente
  % e provocava errore.
  
  % Dunque è da notare che se LAT = 0 potremmo avere due casi:
  % - 0 Reale = LAT = 10log(1);
  % - 0 Finto = per evitare l'errore sul LAT.
  
  if (startattack.pos == 0)  
      LogAttackTime = 0;
      StAtt = 0;
      StPAtt = 0;
  else
  LogAttackTime = log10( ( time_v(stopattack.pos) - time_v(startattack.pos) ) );
  StAtt = time_v(startattack.pos);
  StPAtt = time_v(stopattack.pos);
  end
end
  