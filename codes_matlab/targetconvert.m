function varout = targetconvert(varin)

varin(find(varin==1))=1;
varin(find(varin==2))=2;
varin(find(varin==7))=3;
varin(find(varin==8))=4;
varin(find(varin==3))=3;
varin(find(varin==4))=4;
varin(find(varin==5))=1;
varin(find(varin==6))=2;
varout = varin;

