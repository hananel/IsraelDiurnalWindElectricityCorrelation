
for i=2011:2030
    disp(i);
    [electricityNormalizedYearly,tYearly,electricityNormalizedDaily] = averageIECyear('../YoelCohenData/',1,i);
end