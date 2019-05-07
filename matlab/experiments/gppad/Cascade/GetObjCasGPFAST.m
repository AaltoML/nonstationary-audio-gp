function [Obj,dObjdx] = GetObjCas(x,y,fftCovs,Params,Dims);
  
  % function [Obj,dObjdx] = GetObjCas(x,y,fftCovs,Params,Dims);
  %
  % Computes the transformed envelope derivatives for the
  % cascade version of GPPAD QUICKLY

    [varx,len,mux,varc,vary] = UnpackParamsGP(Params); 
    
    T = Dims.T;
    M = Dims.M;
    Tx = Dims.Tx;

    X = reshape(x,[Tx,M]);  

    %%%%%%%%%%%%%%%%%%%%
    % Likelihood

    A = log(1+exp(X(1:T,:)));
    AComb = prod(A,2);
    
    ObjA = sum(log(A(:))) + sum(1/2*(y./AComb).^2)/varc;
    dAdX = 1./(1+exp(-X(1:T,:)));

    Delta = 1 - 1/varc*y.^2.*AComb.^(-2);
    dObjAdX = [dAdX./A .* (Delta*ones(1,M));zeros(Tx-T,M)];

    %%%%%%%%%%%%%%%%%%%%%
    % Prior
    ObjB = 0;
    dObBdX = zeros(Tx,M);
    
    for m=1:M
      
      xExt = X(:,m)-mux(m);

      fftxExt = fft(xExt);
      
      InvCovX = ifft(fftxExt./fftCovs(:,m));
      
      ObjB = ObjB + 1/(2*varx(m))*InvCovX'*xExt;
      
      dObjBdX(:,m) = InvCovX/varx(m);
    end

    %%%%%%%%%%%%%%%%%%%%%
    % Output
    
    Obj = (ObjA + ObjB);
    dObjdx = (dObjBdX(:)+dObjAdX(:));

    %Obj = ObjA;
    %dObjdx = dObjAdX(:);

    %Obj = ObjB;
    %dObjdx = dObjBdX(:);
