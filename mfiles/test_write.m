tic
fidout = hb_create_file('test.out', 0);
[m,nsta] = size(profiles);
for ii=1:nsta
    err = hb_write_profile(fidout,profiles(ii));
end
disp('finished writing test.out');

fclose(fidout);
clear m nsta ii fidout
toc

