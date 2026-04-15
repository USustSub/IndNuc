function param=convert_thickness(param)
% param.RSb = 2730;             %Rocksalt bottom
param.BZb = param.RSb+param.thickBZ;             %Basal zechstein bottom
param.TBb = param.BZb+param.thickTB;             %Ten Boer bottom
param.SSb = param.TBb+param.thickSS;             %Slochteren sandstone bottom
% param.Cb = 4000;              %Carbonefious bottom
end