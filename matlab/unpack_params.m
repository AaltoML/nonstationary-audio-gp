function [w_sub, w_mod, Wnmf, lik_param] = unpack_params(w,kernel1,kernel2,num_lik_params,D,N,constraints,w_fixed,tune_hypers)

% unwrap the optimisable and fixed parameters with constraints
  w_ind = 1; wf_ind = 1;  
  if tune_hypers(1)
    lik_param = w(1:num_lik_params);
    w_ind = w_ind + num_lik_params;
  else
    lik_param = w_fixed(1:num_lik_params);
    wf_ind = wf_ind + num_lik_params;
  end 
  w_sub = [];
  w_mod = [];
  for i=2:6
    if tune_hypers(i)
      if i<=4
        w_sub = [w_sub; sigmoid(w(w_ind:w_ind+D-1),constraints(i-1,:))];
        w_ind = w_ind + D;
      else
        w_mod = [w_mod; sigmoid(w(w_ind:w_ind+N-1),constraints(i-1,:))];
        w_ind = w_ind + N;
      end
    else
      if i<=4
        w_sub = [w_sub; sigmoid(w_fixed(wf_ind:wf_ind+D-1),constraints(i-1,:))];
        wf_ind = wf_ind + D;
      else
        w_mod = [w_mod; sigmoid(w_fixed(wf_ind:wf_ind+N-1),constraints(i-1,:))];
        wf_ind = wf_ind + N;
      end
    end
  end
  
  if tune_hypers(7)
    Wnmf = reshape(sigmoid(w(w_ind(1):end),constraints(6,:)),[D,N]);
  else
    Wnmf = reshape(sigmoid(w_fixed(wf_ind(1):end),constraints(6,:)),[D,N]);
  end
  
  % calculate conditional variance from marginal
  w_sub(1:D) = w_sub(1:D) .* (1 - w_sub(D+1:2*D));
  w_mod(1:N) = w_mod(1:N) .* (1 - w_mod(N+1:2*N));
  
  % convert lambdas to lengthscales
  w_sub(D+1:2*D) = lambda_map(w_sub(D+1:2*D),kernel1);
  w_mod(N+1:2*N) = lambda_map(w_mod(N+1:2*N),kernel2);

end