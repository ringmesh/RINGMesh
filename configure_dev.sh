#!/bin/sh

if [ ! -f configure.sh ]
then 
    echo ERROR: No configure.sh file.
    exit 1
fi

if [ ! -x configure.sh ]
then 
    echo ERROR: Configure.sh file is not executable.
    exit 1
fi

./configure.sh

echo
echo ======== Code formatter hook =============
if [ ! -x /usr/bin/hg ]
then
    echo ERROR: MercurialI is not instlled on this computer.
    echo If you own this computer install mercurial \(sudo apt-get install mercurial\).
    echo Else ask to your administrator.
    exit 1
fi

if [ ! -d Formatter ]
then
    echo You have no formatter.
    ping -q -c5 192.168.1.4 > /dev/null
    if [[ $? -eq 0 ]]
    then
        echo INFO: Connection to lirac.
        hg clone http://hg/misc/Formatter
    else
        echo ERROR: Impossible to connect to lirac. No Formatter installed.
        exit 1
    fi
fi


grep_result=$(grep hook .hg/hgrc)
if [ -z "${grep_result}" ]
then
    echo INFO: Hook not in hgrc.
    hook_text="\n[hooks]\npre-commit = cd ./Formatter && ./call_uncrustify.sh false ../include/grgmeshlib/ ../src/grgmeshlib/ && cd .."
    hook_text="${hook_text}\n\nchangegroup = cd ./Formatter && hg pull -u && cd .."

    echo ${hook_text} >> .hg/hgrc
    echo INFO: hook added.
fi

if [ ! -x /usr/bin/uncrustify ]
then
    echo WARNING: Uncrustify is not instlled on this computer.
    echo If you own this computer install uncrustify \(sudo apt-get install uncrustify\).
    echo Else ask to your administrator.
fi

if [ ! -x /usr/bin/universalindentgui ]
then
    echo WARNING: UniversalIndentGUI is not instlled on this computer.
    echo If you own this computer install universalindentgui \(sudo apt-get install universalindentgui\).
    echo Else ask to your administrator.
fi

exit 0
