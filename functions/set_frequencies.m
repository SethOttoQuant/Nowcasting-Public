function [frq]  = set_frequencies(string_frq)
  % Set the frequency mix
  % Input:
  %   String frequency: y, q, m, bw, w, d
  %
  % Output:
  %   Number of high frequency periods in each low frequency period to make
  %   the helper matrix J
  k = size(string_frq,1); %number of series
  frq = ones(k,1);
  frq(strcmp('d', string_frq)) = 1;
  frq(strcmp('w', string_frq)) = 7;
  frq(strcmp('bw', string_frq)) = 14;
  frq(strcmp('m', string_frq)) = 28;
  frq(strcmp('q', string_frq)) = 84;
  frq(strcmp('y', string_frq)) = 336;
  min_frq = min(frq);
  frq = frq/min_frq;
  %finalizing
  if min_frq == 1
      frq(frq == 28) = 31;
      frq(frq == 84) = 91;
      frq(frq == 336) = 365;
  elseif min_frq == 7
      frq(frq == 12) = 13;
      frq(frq == 48) = 52;
  elseif min_frq == 14
      frq(frq == 24) = 26;
  end
  return
end
      