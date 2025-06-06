function imaging_scan_extract_bvec_bval(twix)

FullAddr = twix.image.filename;
fileName=strfind(FullAddr,'/');
addr=strcat(FullAddr(1:fileName(end)),'../recon/');
fileName=FullAddr(fileName(end)+1:end-4);

ice=[twix.image.iceParam(7,:)-16384;... % 16384 is offset set by Siemens
    twix.image.iceParam(10:12,:)-16384];
ice(2:4,:)=ice(2:4,:)./sqrt(ice(2,:).^2+ice(3,:).^2+ice(4,:).^2);
ice(isnan(ice))=0;
diff_info=ice(:,1:twix.hdr.Dicom.EchoTrainLength*twix.hdr.Phoenix.sWipMemBlock.alFree{13}:end);
bval=diff_info(1,:);
bvec=diff_info(2:end,:);
write_bvec_bval(bval',bvec',strcat(addr,fileName));
end

function write_bvec_bval(bval,bvec,fileName)

fIDx1=fopen(strcat(fileName,'_bvectors.bvec'),'w');
fIDx2=fopen(strcat(fileName,'_bvalues.bval'),'w');

for i=1:3
    for k=1:size(bvec,1)
        fprintf(fIDx1,"%1.4f ",bvec(k,i));
    end
    fprintf(fIDx1,"\n");
end
for k=1:size(bval,1)
    fprintf(fIDx2,"%d ",bval(k));
end
fclose('all');
end
