 #!/bin/bash
 # Script to launch a push on RINGMeshCoverage repo
 
cd %appveyor_build_folder%
echo [ui] > C:\Users\appveyor\mercurial.ini
echo username = ArnaudBotella >> C:\Users\appveyor\mercurial.ini 
echo [auth] >> C:\Users\appveyor\mercurial.ini
echo bb.prefix = https://bitbucket.org >> C:\Users\appveyor\mercurial.ini
echo bb.username = ArnaudBotella >> C:\Users\appveyor\mercurial.ini
echo bb.password = %password% >> C:\Users\appveyor\mercurial.ini
hg push --new-branch https://bitbucket.org/ring_team/ringmeshcoverage
