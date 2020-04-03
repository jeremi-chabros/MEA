function error_email_fcn(email,errorfilename)

% User input
source = 'mea.matlab.errors@gmail.com';                       %from address (gmail)
destination = email;                                    %to address (any mail service)
myEmailPassword = 'matlaberror123';                           %the password to the 'from' account
subj = 'spike detection failed';         % subject line
msg = {['Bonjour,'],[' '],['This is a message generated automatically from the MATLAB MEA Toolbox spike detection function.'],...
    ['There has been an error in spike detection.'],[' '],...
    [strcat(errorfilename,' — this file has been skipped.')],...
    [' '],['Contact awed2@cam.ac.uk for support if your world is crumbling in around you.']...
    [' '],['Kind regards,']...
    ['The Fish']};            % main body of email.
%set up SMTP service for Gmail
setpref('Internet','E_mail',source);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',source);
setpref('Internet','SMTP_Password',myEmailPassword);
% Gmail server.
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
% Send the email
sendmail(destination,subj,msg);
% [Optional] Remove the preferences (for privacy reasons)
setpref('Internet','E_mail','');
setpref('Internet','SMTP_Server','''');
setpref('Internet','SMTP_Username','');
setpref('Internet','SMTP_Password','');

end