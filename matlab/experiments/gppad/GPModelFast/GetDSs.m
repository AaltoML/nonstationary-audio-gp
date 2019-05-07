function DSs = GetDSs(len,MinLength)

% function DSs = GetDSs(len,MinLength)
  
if len>MinLength
  log2MaxDS = floor(log2(len)-log2(MinLength));  
  MaxDS = 2^log2MaxDS;
  MinDS = 1;
  DSs = round(logspace(log10(MaxDS),log10(MinDS),log2MaxDS+1));
else
  DSs=1;
end
