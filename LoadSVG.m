svgname = 


svgtext=fileread(svgname);
%% Prepare to trasnsform svg string to array
P_status='[MmCcLlSsVvHhZz]';
P_path_str=extractBetween(svgtext,'<path','>');
P_path_str=regexprep(P_path_str,'id="','ID="');
P_d_str=extractBetween(P_path_str,'d="','"');
P_d_str=regexprep(P_d_str,',',' ');
P_d_str=regexprep(P_d_str,'-',' -');
% % P_line_str=regexprep(P_line_str,'(?<!-)[.]',' .');
P_d_str=regexprep(P_d_str,P_status,' $0 ');
% P_line_str=strip(P_line_str);
P_d_str=regexprep(P_d_str,'\s*',' ');
% P_line_str=regexprep(P_line_str,',,',',');