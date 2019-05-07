function OptsNew = ModifyOpts(OptsOld,OptsMod)
  
  % function OptsNew = ModifyOpts(OptsOld,OptsMod)
  %
  % Takes the old options structure and adds / modifies
  % the fields according to OptsMod producing a new
  % options structure OptsNew
  % e.g. if the inputs are:   
  % OptsOld = struct('a',1,'b',[1,2,3],'c','oranges') 
  % OptsMod = struct('b',[1,2,3,5],'c','bananas','d','cheese')
  % OptsNew = struct('a',1,'b',[1,2,3,5],'c','bananas','d','cheese')
    
    OptsNew = OptsOld;
  
    if isstruct(OptsMod)
      FNames = fieldnames(OptsMod);
    
      N = size(FNames);
      
      for n=1:N
        str = ['OptsNew.',FNames{n},' = ','OptsMod.',FNames{n},';'];
        eval(str);
      end
      
    end