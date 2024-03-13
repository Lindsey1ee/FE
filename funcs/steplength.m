function sl = steplength(y_lead, y_trail, GEidx1, GEidx2,Speed, Time, pitch)

% for i = 1:length(Speed)
% sl(i) = ((y_lead(GEidx1(i)) - y_trail(GEidx2(i))) + (Speed(i)*Time(i)))/cosd(pitch) ;
% end

sl = ((y_lead(GEidx1) - y_trail(GEidx2))/cosd(pitch)) + (Speed.*Time) ;
end 

% for i = 1:length(Speed)
% sl = (y_lead(GEidx1) - y_trail(GEidx1)) ;
% end
