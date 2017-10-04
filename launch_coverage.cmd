REM Script to run a push on RINGMeshCoverage repo (and then run code coverage tools)
 
cd %appveyor_build_folder%
echo [ui] > C:\Users\appveyor\mercurial.ini
echo username = ArnaudBotella >> C:\Users\appveyor\mercurial.ini 
echo [auth] >> C:\Users\appveyor\mercurial.ini
echo bb.prefix = https://bitbucket.org >> C:\Users\appveyor\mercurial.ini
echo bb.username = ArnaudBotella >> C:\Users\appveyor\mercurial.ini
echo bb.password = %password% >> C:\Users\appveyor\mercurial.ini
hg push --force --new-branch https://bitbucket.org/ring_team/ringmeshcoverage
