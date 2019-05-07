function [A,X,C,Params,Info] = FBCasGPMAP(Y,Len,Opts) 

  % function [A,X,C,Params,Info] = FBCasGPMAP(Y,Len,Opts) 
  %
  % Cacade demodulation of each column of Y  

[T,D] = size(Y);

[M,D2] = size(Len);

if D2==1
  Len = Len*ones(1,D);
end

A = zeros(T,D,M);
C = zeros(T,D);
X = A;

Params.len = Len;
Params.varc = zeros(1,D);
Params.varx = zeros(M,D);
Params.mux = zeros(M,D);

for d=1:D  
  disp(['Cascade progress: ',num2str(d),'/',num2str(D)])
  lenCur = Len(:,d);
  yCur = Y(:,d);
  
  [ACur,XCur,cCur,ParamsCur,Info] = CasGPMAP(yCur,lenCur, ...
                                             Opts);
  A(:,d,:) = reshape(ACur,[T,1,M]);
  X(:,d,:) = reshape(XCur,[T,1,M]);
  C(:,d) = cCur;

  Params.varc(d) = ParamsCur.varc;
  Params.varx(:,d) = ParamsCur.varx;
  Params.mux(:,d) = ParamsCur.mux;
  
end

A = squeeze(A);
X = squeeze(X);