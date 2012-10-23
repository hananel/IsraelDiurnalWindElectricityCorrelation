im=imread('~/Downloads/North.jpg');
imtool(im);
sz=size(im);

% [203 37 245]
% [196 19 231]
% [183 0 216]
counter = 0;
for i=1:sz(1)
    for j=1:sz(2)
        counter = counter + 1;
        if im(i,j,1)>200 & im(i,j,1)<210 & im(i,j,3)>230 & im(i,j,2)>=30 & im(i,j,2)<40
            pink(i,j)=1;
            I(counter) = i;
            J(counter) = j;
        end
    end
end

imtool(pink);