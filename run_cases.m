cd('G:\My Drive\MP in OSBL\imaging expts\run220613\outputs')
for n = 1:6
    preprocess_MP(n);
end

cd('D:\luci\run220613')
for n = 7:12
    detect_MP_avi(n)
end

cd('G:\My Drive\MP in OSBL\imaging expts\run220613\outputs')
for n = 7:12
    preprocess_MP(n);
end

cd('D:\luci\run220613')
for n = 13:18
    detect_MP_avi(n)
end

cd('G:\My Drive\MP in OSBL\imaging expts\run220613\outputs')
for n = 13:18
    preprocess_MP(n);
end