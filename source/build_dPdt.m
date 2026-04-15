function [dPdt]=build_dPdt(Ny,z,t2,param)
    if (t2<0)
        dPdt.L=zeros(Ny,1);
        dPdt.R=zeros(Ny,1);
    else
    TBt=param.BZb-param.ytop;
    SSt=param.TBb-param.ytop;
    SSb=param.SSb-param.ytop;
    offset=param.offset;
    
    dPdtSS=param.dPdtSS;
    dPdtTB=param.dPdtTB;
    Ddiff=param.Ddiff;
    hdiff=sqrt(Ddiff*t2);
    
    dPdt.L=nan(Ny,1);
    if (dPdtTB<0)
    dPdt.L(z<=TBt-offset-hdiff,1)=0;
    dPdt.L(z>TBt-offset+hdiff&z<=SSt-offset-hdiff,1)=dPdtTB;
    else
        dPdt.L(z<=SSt-offset-hdiff,1)=0;
    end
    dPdt.L(z>SSt-offset+hdiff&z<=SSb-offset-hdiff,1)=dPdtSS;
    dPdt.L(z>SSb-offset+hdiff,1)=0;
    dPdt.L=fillmissing(dPdt.L,'linear');
    
%     d2Pdt.L=[diff(dPdt.L);0];
    
    dPdt.R=nan(Ny,1);
    if (dPdtTB<0)
    dPdt.R(z<=TBt-hdiff,1)=0;
    dPdt.R(z>TBt+hdiff&z<=SSt-hdiff,1)=dPdtTB;
    else
        dPdt.R(z<=SSt-hdiff,1)=0;
    end
    dPdt.R(z>SSt+hdiff&z<=SSb-hdiff,1)=dPdtSS;
    dPdt.R(z>SSb+hdiff,1)=0;
    dPdt.R=fillmissing(dPdt.R,'linear');
    
%     d2Pdt.R=[diff(dPdt.R);0];
    end
end

