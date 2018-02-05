#!/bin/bash
# create the list for program shapestring, 

# ChangeLog 2008-04-28
#    add the option --pdbaa
#    add the option --pdbaatype
# ChangeLog 2009-12-21
#    redirect the error message to stderr with >&2
#    
# Copyright© Nanjiang Shu 
# Department of Materials and Environmental Chemistry, Stockholm University,
# Sweden
# Email: nanjiang.shu@mmk.su.se

usage="
Usage: cShapeStringList.sh [options] idListFile
Create the list file for the program 'shapestring'
options:
  -n                  : do not check the existence of necessary programs
  -pdbaa|--pdbaa path : set the pdbaa path, default = $DATADIR/pdbaa
  -pdbaatype 0|1      : 0 -- e.g. 1bmf_A.aa  1 -- e.g. 1BMFA.aa, default = 1
  -k|--keep           : do not delete temporary files
  -h|--help           : print this help message and exit

Created 2008-02-04, updated 2010-05-27, Nanjiang Shu
"

function PrintHelp()
{
    echo "$usage"
}

function IsProgExist()#{{{
# usage: IsProgExist prog
# prog can be both with or without absolute path
{
    type -P $1 &>/dev/null || { echo "The program \"$1\" is required but it's not installed. Aborting." >&2; exit 1; }
}
#}}}
if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

isCheckExistence=true
isKeepTempFile=false
pdbaapath=
pdbaatype=1 #changed 2009-12-21
binpath=$BUILDSEQMAP/bin

while [ "$1" != "" ]; do
    case $1 in
        -h|--help)PrintHelp;exit;;
        -n )isCheckExistence=false;;
        -pdbaa|--pdbaa)pdbaapath=$2;shift;;
        -pdbaatype |--pdbaatype)pdbaatype=$2;shift;;
        -k|--keep )isKeepTempFile=true;;
        *) idListFile=$1;;
    esac
    shift
done

if [ ! -f $idListFile ]; then
    echo "Error! idListFile = \"$idListFile\" can not be opend" >&2
    exit
fi

# check the necessary programs
if [ "$isCheckExistence" == "true" ]; then
    IsProgExist $binpath/getpdbfilepath
    IsProgExist $binpath/getpdbaafilepath
    IsProgExist $binpath/getdsspfilepath
    IsProgExist awk
    IsProgExist paste
fi

tmpPDBIDListFile=$idListFile.$$.pdblistfile
tmpPDBAAListFile=$idListFile.$$.pdbaalistfile
tmpDSSPListFile=$idListFile.$$.dssplistfile
tmpChainIDListFile=$idListFile.$$.chainIDList
errFile=$idListFile.$$.err

if [ -f "$idListFile" ]; then 
    $binpath/getpdbfilepath -l $idListFile 1> $tmpPDBIDListFile 2> $errFile
    if [ -s $errFile ]; then 
        cat $errFile >&2
        exit
    fi
    if [ "$pdbaatype" ==  "0" ] ; then
        $binpath/getpdbaafilepath -d $pdbaapath -l $idListFile 1> $tmpPDBAAListFile 2> $errFile
        if [ -s $errFile ]; then 
            cat $errFile  >&2
            exit
        fi
    else
        awk -v path=$pdbaapath '{print path "/" $1 ".aa"}' $idListFile > $tmpPDBAAListFile
    fi

    $binpath/getdsspfilepath -l $idListFile 1> $tmpDSSPListFile 2> $errFile
    if [ -s $errFile ]; then 
        cat $errFile >&2
        exit
    fi
    awk '{print "C"substr($1,5,1)}' $idListFile > $tmpChainIDListFile 

    paste $tmpPDBIDListFile $tmpChainIDListFile $idListFile $tmpPDBAAListFile $tmpDSSPListFile

    if [ "$isKeepTempFile" == "false" ]; then
        rm -f $tmpPDBIDListFile
        rm -f $tmpChainIDListFile
        rm -f $tmpPDBAAListFile
        rm -f $tmpDSSPListFile
        rm -f $errFile
    fi
fi

