function [lout] = lambda_map(lin, kernel)

  switch kernel
    case 'exp'
      lout = 1 ./ lin;
    case 'matern32'
      lout = sqrt(3) ./ lin;
    case 'matern52'
      lout = sqrt(5) ./ lin;
    case 'matern72'
      lout = sqrt(7) ./ lin;
    case 'se'
      error('what should the lambda mapping be for se? - not sure')
  end
