#!/bin/bash

echo Cpp fomatting is started.

# If no parameter, error
if [[ $# -lt 3 ]]
then
    # TODO a updater
    echo "Syntax is: call_uncrustify.sh uncrustify.cfg create_backup file_or_dir1 file_or_dir2 ..."
    exit 1
fi

cpp_extensions=( "cpp" "h" "hpp" )
#faudrait virer les moc...
# en fait en param faudrait prendre une option qui pourrait etre une expression reguliere et on formatte pas les fichiers qui correspondent a ce format
for arg in "$@"
do
    if [[ ${arg} == ${0} ]]
    then
        continue
    elif [[ ${arg} == ${1} ]]
    then
        if [[ ! -f ${1}  ]]
        then
            echo ERROR uncrustify config file does not exist.
            exit 1
        else
            continue
        fi
    elif [[ ${arg} == ${2} ]]
    then
        if [[ ${arg} != "true" && ${arg} != "false" ]]
        then
            echo ERROR second parameter is not a boolean.
            exit 1
        fi
        continue
    fi

    if [[ -d ${arg} ]]
    then
        for filesuffix in "${cpp_extensions[@]}"
        do
            file_list=`find ${arg} -name "*.${filesuffix}" -type f`
            for file2indent in $file_list
            do 
                uncrustify -q -f "$file2indent" -c ${1} -o ${file2indent}.tmp &> /dev/null
                if [[ -f ${file2indent}.tmp ]]
                then
                    if [[ ${2} == "true" ]]
                    then
                        mv "$file2indent" "${file2indent}.bak"
                    fi
                    diff_exists=$(diff ${file2indent}.tmp $file2indent)
                    if [[ ${diff_exists} != "" ]]
                    then
                        echo File "$file2indent" has been formatted.
                    fi
                    mv ${file2indent}.tmp "$file2indent"
                fi
            done
        done
    elif [[ -f "${arg}" ]]
    then
        echo "Indenting one file ${arg}"
        uncrustify -q -f "$arg" -c ${1} -o ${arg}.tmp
        if [[ -f ${file2indent}.tmp ]]
        then
            if [[ ${2} == "true" ]]
            then
                mv "$arg" "${arg}.bak"
            fi
            diff_exists=$(diff ${file2indent}.tmp $file2indent)
            if [[ ${diff_exists} != "" ]]
            then
                echo File "$file2indent" has been formatted.
            fi

            mv ${arg}.tmp "$arg"
        fi
    else
        echo ERROR input file/directory does not exist.
        exit 1
    fi
done

echo Cpp formatting is finished.
exit 0
