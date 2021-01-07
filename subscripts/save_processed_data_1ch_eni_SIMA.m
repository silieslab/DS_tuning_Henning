function save_processed_data_1ch_eni_SIMA(in,file);



strct=in;
if(isfield(in,'ch1b'))
    strct = rmfield(strct,'ch1b');
    strct = rmfield(strct,'ch3');

end

pathroot=file;

f=find((pathroot == '\') + (pathroot == '/'));
f=f(end);
pathroot = pathroot((f+1):end); 
flyroot = in.dataID;
save([flyroot '_' pathroot '_pData_SIMA_only_m.mat'],'strct', '-v7.3' );

end